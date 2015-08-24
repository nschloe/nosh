// @HEADER
//
//    Builds the kinetic energy operator and its variants.
//    Copyright (C) 2010--2012  Nico Schlömer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
// =============================================================================
// includes
#include "parameter_matrix_dkeo_dp.hpp"

#include <map>
#include <string>

#include "mesh.hpp"
#include "scalar_field_base.hpp"
#include "vector_field_base.hpp"

#include <Tpetra_Vector.hpp>
#include <Tpetra_Import.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Teuchos_Tuple.hpp>

#include <ml_epetra_preconditioner.h>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif

namespace nosh
{
namespace parameter_matrix
{
// =============================================================================
DkeoDP::
DkeoDP(
    const std::shared_ptr<const nosh::mesh> &mesh,
    const std::shared_ptr<const nosh::scalar_field::base> &thickness,
    const std::shared_ptr<nosh::vector_field::base> &mvp,
    const std::string & param_name
   ):
  parameter_object(),
  Tpetra::CrsMatrix<double,int,int>(mesh->build_complex_graph()),
  mesh_(mesh),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  keo_fill_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: DkeoDP::refill_")),
#endif
  thickness_(thickness),
  mvp_(mvp),
  alpha_cache_(),
  alpha_cache_up_to_date_(false),
  param_name_(param_name)
{
}
// =============================================================================
DkeoDP::
~DkeoDP()
{
}
// =============================================================================
const std::map<std::string, double>
DkeoDP::
get_parameters() const
{
  return mvp_->get_parameters();
}
// =============================================================================
void
DkeoDP::
refill_(const std::map<std::string, double> & params)
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*keo_fill_time_);
#endif
  this->resumeFill();

  mvp_->set_parameters(params);

  this->setAllToScalar(0.0);

#ifndef NDEBUG
  TEUCHOS_ASSERT(mesh_);
#endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
#ifndef NDEBUG
  TEUCHOS_ASSERT(thickness_);
  TEUCHOS_ASSERT(mvp_);
#endif

  const std::vector<edge> edges = mesh_->my_edges();
  if (!alpha_cache_up_to_date_) {
    this->build_alpha_cache_(edges, mesh_->edge_coefficients());
  }

  std::vector<double> v(3);

  // Loop over all edges.
  for (std::size_t k = 0; k < edges.size(); k++) {
    // Compute the integral
    //
    //    I = \int_{x0}^{xj} (xj-x0).A(x) / |xj-x0| dx
    //
    // numerically by the midpoint rule, i.e.,
    //
    //    I ~ |xj-x0| * (xj-x0) . A(0.5*(xj+x0)) / |xj-x0|.
    //
    // Project vector field onto the edge.
    // Instead of first computing the projection over the normalized edge and
    // then multiply it with the edge length, don't normalize the edge vector.
    double a_int = mvp_->get_edge_projection(k);
    // param_name_ is set in the KEO building routine.
    double dAdPInt = mvp_->get_d_edge_projection_dp(k, param_name_);
    double sin_a_int, cos_a_int;
    sincos(a_int, &sin_a_int, &cos_a_int);
    v[0] =  dAdPInt * sin_a_int;
    v[1] = -dAdPInt * cos_a_int;
    v[2] = 0.0;
    // We'd like to insert the 2x2 matrix
    //
    //     [   alpha                   , - alpha * exp(-IM * a_int) ]
    //     [ - alpha * exp(IM * a_int),   alpha                       ]
    //
    // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
    // that shares and edge.
    // Do that now, just blockwise for real and imaginary part.
    v[0] *= alpha_cache_[k];
    v[1] *= alpha_cache_[k];
    v[2] *= alpha_cache_[k];
    auto vals = Teuchos::tuple(
      Teuchos::tuple(v[2],  0.0,   v[0], v[1]),
      Teuchos::tuple( 0.0,  v[2], -v[1], v[0]),
      Teuchos::tuple(v[0], -v[1],  v[2],  0.0),
      Teuchos::tuple(v[1],  v[0],   0.0, v[2])
      );
    const Teuchos::Tuple<int,4> & idx = mesh_->edge_gids_complex[k];
    for (int i = 0; i < 4; i++) {
      int num = this->sumIntoGlobalValues(idx[i], idx, vals[i]);
#ifndef NDEBUG
      TEUCHOS_ASSERT_EQUALITY(num, 4);
#endif
    }
  }

  this->fillComplete();

  return;
}
// =============================================================================
void
DkeoDP::
build_alpha_cache_(
    const std::vector<edge> & edges,
    const std::vector<double> & edge_coefficients
    ) const
{
  // This routine serves the one and only purpose of caching the thickness
  // average. The cache is used in every call to this->fill().  This is
  // somewhat problematic since the map of V is principally not known here.
  // Also, it is typically a nonoverlapping map whereas some edges do sit on a
  // processor boundary, so actually the values of V are needed in an
  // overlapping map.
  // Fair enough. Let's distribute the vales of V to an overlapping map here.
  alpha_cache_ = std::vector<double>(edges.size());

  std::map<std::string, double> dummy;
  const auto thickness_values = thickness_->get_v(dummy);

  auto overlapMap = mesh_->overlap_map();
  // We need to make sure that thickness_values are distributed on the overlap
  // map.
  // Make sure to use Import here instead of Export as the vector that we want
  // to build is overlapping, "larger". If the "smaller", non-overlapping
  // vector is exported, only the values on the overlap would only be set on
  // one processor.
  Tpetra::Vector<double,int,int> thicknessOverlap(Teuchos::rcp(overlapMap));
  auto importer = Tpetra::createImport(
      thickness_values.getMap(),
      Teuchos::rcp(overlapMap)
      );
  thicknessOverlap.doImport(thickness_values, *importer, Tpetra::INSERT);

  auto t_data = thicknessOverlap.getData();

  int gid0, gid1;
  int lid0, lid1;
  for (std::size_t k = 0; k < edges.size(); k++) {
    // Get the ID of the edge endpoints in the map of get_v(). Well...
    gid0 = mesh_->gid(std::get<0>(edges[k]));
    lid0 = overlapMap->getLocalElement(gid0);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid0 < 0,
        "The global index " << gid0
        << " does not seem to be present on this node."
       );
#endif
    gid1 = mesh_->gid(std::get<1>(edges[k]));
    lid1 = overlapMap->getLocalElement(gid1);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid1 < 0,
        "The global index " << gid1
        << " does not seem to be present on this node."
       );
#endif
    // Update cache.
    alpha_cache_[k] = edge_coefficients[k]
      * 0.5 * (t_data[lid0] + t_data[lid1]);
  }

  alpha_cache_up_to_date_ = true;

  return;
}
// =============================================================================
}  // namespace parameter_matrix
}  // namespace nosh
