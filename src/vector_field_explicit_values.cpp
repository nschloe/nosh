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

  // TODO(nschloe): make mb_ private again
  moab::Range verts = mesh.mbw_->get_entities_by_dimension(0, 0);

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
    // values at the adjacent vertices.
    Eigen::Vector3d av = 0.5 * (A0 + A1);

    const auto edge_coords = mesh.get_coords(v0) - mesh.get_coords(v1);

    edgeProjectionCache_[k] = av.dot(edge_coords);
  }
  return;
}
// ============================================================================
explicit_values::
~explicit_values()
= default;
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
get_scalar_parameters() const
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
  if (dParamName.compare("mu") == 0) {
    return edgeProjectionCache_[edge_index];
  } else {
    return 0.0;
  }
}
// ============================================================================
} // namespace vector_field
} // namespace nosh
