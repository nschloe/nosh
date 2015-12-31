// @HEADER
//
//    Helper class for writing statistics and states.
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

#include "observer.hpp"

#include <string>

#include "model_evaluator_base.hpp"
#include "mesh.hpp"

namespace nosh
{
// ============================================================================
observer::
observer(
    const std::shared_ptr<const nosh::model_evaluator::base> &model_eval,
    const std::string & csv_filename,
    const std::string & cont_param_name,
    const bool is_turning_point_continuation
    ) :
  model_eval_(model_eval),
  csv_writer_(csv_filename, " "),
  cont_param_name_(cont_param_name),
  is_turning_point_continuation_(is_turning_point_continuation)
{
}
// ============================================================================
observer::
~observer ()
{
}
// ============================================================================
void
observer::
observeSolution(const Thyra::VectorBase<double> &soln)
{
  // TODO
  //model_eval_->mesh()->insert(soln, "psi");
  //model_eval_->mesh()->write(0.0);

  return;
}
// ============================================================================
void
observer::
observeSolution(
  const Thyra::VectorBase<double> & soln,
  double param_val
  )
{
  // This if-else hack is necessary as different continuation algorithms
  // call printSolution() a different number of times per step, e.g.,
  // to store solutions, null vectors, and so forth.
  if (is_turning_point_continuation_) {
    this->observe_turning_point_continuation_(soln, param_val);
  } else {
    this->observeContinuation_(soln, param_val);
  }

  return;
}
// ============================================================================
void
observer::
observeContinuation_(
    const Thyra::VectorBase<double>  &soln,
    const double param_val
    )
{
  static int index = -1;
  index++;

  this->save_continuation_statistics_(soln, param_val, index);

  // Storing the parameter value as "time" variable here is convenient, but
  // has a downside: The default output format ExodusII insists that the
  // values for time are monotonically increasing. The parameter, however,
  // can decrease. Since there's no hard reason for the monotonicity condition,
  // many things will continue to work fine if the time data isn't monotonous.
  // The display in ParaView is one example where it doesn't work so well.
  // As a work-around for that, param_val could be replaced by index.
  // TODO
  //model_eval_->mesh()->insert(soln, "psi");
  //model_eval_->mesh()->write(param_val);

  return;
}
// ============================================================================
void
observer::
observe_turning_point_continuation_(
    const Thyra::VectorBase<double> &soln,
    const double param_val
    )
{
  static int index = -1;
  static bool is_solution = false;

  // alternate between solution and nullvector
  is_solution = !is_solution;
  if (is_solution) {
    index++;
    this->save_continuation_statistics_(soln, param_val, index);
    // TODO
    //model_eval_->mesh()->insert(soln, "psi");
    //model_eval_->mesh()->write(index);
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not yet implemented.");
  }
  // This part of the code used to write state and null vector alternately
  // for turning point continuation, but because of how Mesh is
  // organized, it seems impossible to first write to one file, then to
  // another with with the same mesh. Need to investigate.
  return;
}
// ============================================================================
void
observer::
save_continuation_statistics_(
    const Thyra::VectorBase<double> &soln,
    const double param_val,
    const int step_index
    )
{
  // Construct parameter list to stuff into the csv_writer_.
  Teuchos::ParameterList paramList;
  paramList.set("(0) step", step_index);

  // Continuation parameter.
  paramList.set("(1) "+cont_param_name_, param_val);

  // Some extra stats.
  paramList.set("(2) Gibbs energy", model_eval_->gibbs_energy(soln));
  paramList.set("(2) ||x||_2 scaled", model_eval_->norm(soln));

  // Write out header.
  if (step_index == 0)
    csv_writer_.write_header(paramList);
  // Write out the data.
  csv_writer_.write_row(paramList);
  return;
}
// ============================================================================
} // namespace nosh
