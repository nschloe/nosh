// @HEADER
//
//    Builds the kinetic energy operator and its variants.
//    Copyright (C) 2010--2012  Nico Schl\"omer
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
#include "nosh/ParameterMatrix_Keo.hpp"

#include <map>
#include <string>

#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Virtual.hpp"
#include "nosh/VectorField_Virtual.hpp"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>

#include <ml_epetra_preconditioner.h>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif

namespace Nosh
{
namespace ParameterMatrix
{
// =============================================================================
Keo::
Keo(const Teuchos::RCP<const Nosh::StkMesh> &mesh,
    const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness,
    const Teuchos::RCP<Nosh::VectorField::Virtual> &mvp
   ):
  Virtual(mesh),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  keoFillTime_(Teuchos::TimeMonitor::getNewTimer("Nosh: Keo::fill_")),
#endif
  thickness_(thickness),
  mvp_(mvp),
  alphaCache_(),
  alphaCacheUpToDate_(false)
{
}
// =============================================================================
Keo::
~Keo()
{
}
// =============================================================================
Teuchos::RCP<Virtual>
Keo::
clone() const
{
  return Teuchos::RCP<Keo>(new Keo(*this));
}
// =============================================================================
const std::map<std::string, double>
Keo::
getParameters() const
{
  return mvp_->getParameters();
}
// =============================================================================
void
Keo::
refill_(const std::map<std::string, double> & params)
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*keoFillTime_);
#endif

  mvp_->setParameters(params);

  // Zero-out the matrix.
  TEUCHOS_ASSERT_EQUALITY(0, this->PutScalar(0.0));

#ifndef NDEBUG
  TEUCHOS_ASSERT(!mesh_.is_null());
#endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
#ifndef NDEBUG
  TEUCHOS_ASSERT(!thickness_.is_null());
  TEUCHOS_ASSERT(!mvp_.is_null());
#endif

  const Teuchos::Array<edge> edges = mesh_->getEdgeNodes();
  if (!alphaCacheUpToDate_)
    this->buildAlphaCache_(edges, mesh_->getEdgeCoefficients());

  double v[3];
  Epetra_SerialDenseMatrix A(4, 4);
  // Loop over all edges.
  for (auto k = 0; k < edges.size(); k++) {
    // ---------------------------------------------------------------
    // Compute the integral
    //
    //    I = \int_{x0}^{xj} (xj-x0).A(x) / |xj-x0| dx
    //
    // numerically by the midpoint rule, i.e.,
    //
    //    I ~ |xj-x0| * (xj-x0) . A(0.5*(xj+x0)) / |xj-x0|.
    //
    // -------------------------------------------------------------------
    // Project vector field onto the edge.
    // Instead of first computing the projection over the normalized edge
    // and then multiply it with the edge length, don't normalize the
    // edge vector.
    // Fill v with
    // Re(-exp(i Aint))
    // Im(-exp(i Aint))
    // 1.0
    const double aInt = mvp_->getEdgeProjection(k);
    // If compiled with GNU (and maybe other compilers), we could use
    // sincos() here to compute sin and cos simultaneously.
    // PGI, for one, doesn't support sincos, though.
    v[0] = -cos(aInt);
    v[1] = -sin(aInt);
    v[2] = 1.0;
    //sincos(aInt, &sinAInt, &cosAInt);
    // We'd like to insert the 2x2 matrix
    //
    //     [   alpha                   , - alpha * exp(-IM * aInt) ]
    //     [ - alpha * exp(IM * aInt),   alpha                       ]
    //
    // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
    // that shares and edge.
    // Do that now, just blockwise for real and imaginary part.
    v[0] *= alphaCache_[k];
    v[1] *= alphaCache_[k];
    v[2] *= alphaCache_[k];
    A(0, 0) = v[2];
    A(0, 1) = 0.0;
    A(0, 2) = v[0];
    A(0, 3) = v[1];
    A(1, 0) = 0.0;
    A(1, 1) = v[2];
    A(1, 2) = -v[1];
    A(1, 3) = v[0];
    A(2, 0) = v[0];
    A(2, 1) = -v[1];
    A(2, 2) = v[2];
    A(2, 3) = 0.0;
    A(3, 0) = v[1];
    A(3, 1) = v[0];
    A(3, 2) = 0.0;
    A(3, 3) = v[2];
    TEUCHOS_ASSERT_EQUALITY(
        0,
        this->SumIntoGlobalValues(mesh_->globalIndexCache[k], A)
        );
    // -------------------------------------------------------------------
  }

  // calls FillComplete by default
  TEUCHOS_ASSERT_EQUALITY(0, this->GlobalAssemble());
  return;
}
// =============================================================================
void
Keo::
buildAlphaCache_(
    const Teuchos::Array<edge> & edges,
    const Teuchos::ArrayRCP<const double> &edgeCoefficients
    ) const
{
  // This routine serves the one and only purpose of caching the
  // thickness average. The cache is used in every call to this->fill().
  // This is somewhat problematic since the map of V is principally
  // not known here. Also, it is typically a nonoverlapping map whereas
  // some edges do sit on a processor boundary, so actually the values
  // of V are needed in an overlapping map.
  // Fair enough. Let's distribute the vales of V to an overlapping
  // map here.
  alphaCache_ = Teuchos::ArrayRCP<double>(edges.size());

  std::map<std::string, double> dummy;
  const Epetra_Vector thicknessValues = thickness_->getV(dummy);

  Teuchos::RCP<const Epetra_Map> overlapMap = mesh_->getNodesOverlapMap();
  // We need to make sure that thicknessValues are distributed on
  // the overlap map.
  // Make sure to use Import here instead of Export as the vector
  // that we want to build is overlapping, "larger". If the "smaller",
  // non-overlapping vector is exported, only the values on the overlap
  // would only be set on one processor.
  Epetra_Vector thicknessOverlap(*overlapMap);
  Epetra_Import importer(*overlapMap, thicknessValues.Map());
  TEUCHOS_ASSERT_EQUALITY(
      0,
      thicknessOverlap.Import(thicknessValues, importer, Insert)
      );

  Teuchos::Tuple<int, 2> gid;
  Teuchos::Tuple<int, 2> lid;
  for (Teuchos::Array<edge>::size_type k = 0;
       k < edges.size();
       k++) {
    // Get the ID of the edge endpoints in the map of
    // getV(). Well...
    gid[0] = mesh_->gid(std::get<0>(edges[k]));
    lid[0] = overlapMap->LID(gid[0]);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid[0] < 0,
        "The global index " << gid[0]
        << " does not seem to be present on this node."
       );
#endif
    gid[1] = mesh_->gid(std::get<1>(edges[k]));
    lid[1] = overlapMap->LID(gid[1]);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid[1] < 0,
        "The global index " << gid[1]
        << " does not seem to be present on this node."
       );
#endif
    // Update cache.
    alphaCache_[k] = edgeCoefficients[k]
      * 0.5 * (thicknessOverlap[lid[0]] + thicknessOverlap[lid[1]]);
  }

  alphaCacheUpToDate_ = true;

  return;
}
// =============================================================================
}  // namespace ParameterMatrix
}  // namespace Nosh
