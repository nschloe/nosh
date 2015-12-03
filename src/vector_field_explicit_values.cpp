// @HEADER
//
//    Query routines for a vector potential with explicitly given values.
//    Copyright (C) 2012  Nico Schlömer
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

#include "vector_field_explicit_values.hpp"

#include <map>
#include <string>

#include "mesh.hpp"

namespace nosh
{
namespace vector_field
{
// ============================================================================
explicit_values::
explicit_values(
    const nosh::mesh & mesh,
    const std::string & field_name,
    const double mu
    ) :
  mu_(mu),
  edgeProjectionCache_(mesh.my_edges().size())
{
  // Initialize the cache.
  const std::vector<edge> edges = mesh.my_edges();

  moab::ErrorCode ierr;
  moab::Range verts;
  // TODO make mb_ private again
  ierr = mesh.mb_->get_entities_by_dimension(0, 0, verts);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  const auto data = mesh.get_data(field_name, verts);

  // Loop over all edges and create the cache.
  for (std::size_t k = 0; k < edges.size(); k++) {
    const auto v0 = std::get<0>(edges[k]);
    const auto v1 = std::get<1>(edges[k]);

    const auto idx0 = mesh.local_index(v0);
    const auto idx1 = mesh.local_index(v1);

    Eigen::Vector3d A0(data[3*idx0], data[3*idx0 + 1], data[3*idx0 + 2]);
    Eigen::Vector3d A1(data[3*idx1], data[3*idx1 + 1], data[3*idx1 + 2]);

    // Approximate the value at the midpoint of the edge by the average of the
    // values at the adjacent nodes.
    Eigen::Vector3d av = 0.5 * (A0 + A1);

    std::vector<double> coords0 = mesh.get_coords(v0);
    std::vector<double> coords1 = mesh.get_coords(v1);

    Eigen::Vector3d edge_coords(
        coords0[0] - coords1[0],
        coords0[1] - coords1[1],
        coords0[2] - coords1[2]
        );

    edgeProjectionCache_[k] = av.dot(edge_coords);
  }

// TODO resurrect this
//#ifndef NDEBUG
#if 0
  // Do a quick sanity check for the edgeProjectionCache_.  It happens too
  // often that the reader elements aren't specified correctly and stk_io
  // *silently* "reads" only zeros.  Use the fake logical "isNonzeroLocal"
  // since Teuchos::Comm<int> doesn't have logical any() or all() operations.
  bool isZeroLocal = true;
  for (std::size_t k = 0; k < edgeProjectionCache_.size(); k++) {
    if (fabs(edgeProjectionCache_[k]) > 1.0e-10) {
      isZeroLocal = false;
      break;
    }
  }
  bool isZeroGlobal;
  Teuchos::reduceAll(
      *mesh.getComm(),
      Teuchos::REDUCE_AND,
      1,
      &isZeroLocal,
      &isZeroGlobal
      );

  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      isZeroGlobal,
      "Field \"" << field_name << "\" seems empty. Was it read correctly?"
      );
#endif
  return;
}
// ============================================================================
explicit_values::
~explicit_values()
{
}
// ============================================================================
void
explicit_values::
set_parameters(const std::map<std::string, double> & params)
{
  mu_ = params.at("mu");
  return;
}
// ============================================================================
const std::map<std::string, double>
explicit_values::
get_parameters() const
{
  return {{"mu", mu_}};
}
// ============================================================================
double
explicit_values::
get_edge_projection(const unsigned int edge_index) const
{
  return mu_ * edgeProjectionCache_[edge_index];
}
// ============================================================================
double
explicit_values::
get_d_edge_projection_dp(
    const unsigned int edge_index,
    const std::string & dParamName
    ) const
{
  if (dParamName.compare("mu") == 0)
    return edgeProjectionCache_[edge_index];
  else
    return 0.0;
}
// ============================================================================
} // namespace vector_field
} // namespace nosh
