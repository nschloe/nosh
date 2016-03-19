#include "parameter_object.hpp"

namespace nosh
{
// ============================================================================
parameter_object::
parameter_object():
  build_parameters_()
{
}
// ============================================================================
parameter_object::
~parameter_object()
{
}
// ============================================================================
void
parameter_object::
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
} // namespace nosh
