#include "parameter_object.hpp"

namespace nosh
{
// ============================================================================
void
parameter_object::
set_parameters(
    const std::map<std::string, double> & scalar_params,
    const std::map<std::string, std::shared_ptr<const Tpetra::Vector<double, int, int>>> & vector_params
    )
{
  // Cache the construction of the matrix.
  // This is useful because in the continuation context, the matrix is called a
  // number of times with the same arguments (in compute_f, getJacobian(), and
  // get_preconditioner().
  bool needs_refill;
  if (build_parameters_scalar_.empty() || build_parameters_vector_.empty()) {
    needs_refill = true;
  }
  else if (!vector_params.empty()) {
    // don't bother checking vector params for now. always flag for rebuild.
    // check vector parameters
    needs_refill = true;
  } else {
    needs_refill = false;
    for (auto const &param: build_parameters_scalar_) {
      if (param.second != scalar_params.at(param.first)) {
        needs_refill = true;
        break;
      }
    }
  }

  if (needs_refill) {
    this->refill_(scalar_params, vector_params);
    // set build_parameters_*_
    for (auto &param: build_parameters_scalar_) {
      param.second = scalar_params.at(param.first);
    }
    for (auto &param: build_parameters_vector_) {
      param.second = vector_params.at(param.first);
    }
  }

  return;
}
// ============================================================================
} // namespace nosh
