#ifndef NOSH_NONLINEAR_SOLVER_HPP
#define NOSH_NONLINEAR_SOLVER_HPP

#include <Teuchos_ParameterList.hpp>
#include <Thyra_ModelEvaluatorDefaultBase.hpp>
#include <Tpetra_Vector.hpp>

#include <boost/any.hpp>

#include <map>
#include <memory>

using list = std::map<std::string, boost::any>;
namespace nosh {
  void
  nonlinear_solve(
      const std::shared_ptr<Thyra::ModelEvaluatorDefaultBase<double>> & model,
      std::map<std::string, boost::any> solver_params
      );
}
#endif // NOSH_NONLINEAR_SOLVER_HPP
