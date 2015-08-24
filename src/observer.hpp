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
#ifndef NOSH_NOXOBSERVER_H
#define NOSH_NOXOBSERVER_H
// =============================================================================
#include <string>

#include <Teuchos_RCP.hpp>
#include <Piro_ObserverBase.hpp>

#include "csv_writer.hpp"
// =============================================================================
// forward declarations
namespace nosh
{
namespace model_evaluator
{
class base;
}
}
// =============================================================================
namespace nosh
{

class observer: public Piro::ObserverBase<double>
{
public:
  //! Constructor
  observer(
      const std::shared_ptr<const nosh::model_evaluator::base> &model_eval,
      const std::string & csv_filename = "",
      const std::string & cont_param_name = "",
      const bool is_turning_point_continuation = false
      );

  //! Destructor
  virtual
  ~observer ();

  virtual
  void
  observeSolution(const Thyra::VectorBase<double> &soln);

  virtual
  void
  observeSolution(
      const Thyra::VectorBase<double> & soln,
      double param_val
      );

protected:
private:
  void
  observeContinuation_(
      const Thyra::VectorBase<double> &soln,
      const double param_val
      );

  void
  observe_turning_point_continuation_(
      const Thyra::VectorBase<double> &soln,
      const double param_val
      );

  void
  save_continuation_statistics_(
      const Thyra::VectorBase<double> &soln,
      const double param_val,
      const int step_index
      );

private:
  const std::shared_ptr<const nosh::model_evaluator::base> model_eval_;
  nosh::csv_writer csv_writer_;
  const std::string cont_param_name_;
  const bool is_turning_point_continuation_;
};
} // namespace nosh
#endif // NOSH_NOXOBSERVER_H
