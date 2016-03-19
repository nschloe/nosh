#include "parameter_matrix_base.hpp"

#include "mesh.hpp"

namespace nosh
{
namespace parameter_matrix
{
// ============================================================================
base::
base(const std::shared_ptr<const nosh::mesh> &mesh):
  Tpetra::CrsMatrix<double,int,int>(mesh->build_complex_graph()),
  mesh_(mesh),
  build_parameters_()
{
}
// ============================================================================
base::
~base()
{
}
// ============================================================================
void
base::
set_parameters(const std::map<std::string, double> &params)
{
  // Cache the construction of the matrix.
  // This is useful because in the continuation context, the matrix is called a
  // number of times with the same arguments (in compute_f, getJacobian(), and
  // get_preconditioner().
  bool needs_refill;
  if (build_parameters_.empty()) {
    needs_refill = true;
  } else {
    needs_refill = false;
    for (auto const &build_param: build_parameters_) {
      if (build_param.second != params.at(build_param.first)) {
        needs_refill = true;
        break;
      }
    }
  }

  if (needs_refill) {
    this->refill_(params);
    for (auto &build_param: build_parameters_) {
      build_param.second = params.at(build_param.first);
    }
  }

  return;
}
// ============================================================================
} // namespace parameter_matrix
} // namespace nosh
