#ifndef NOSH_NONLINEAR_SOLVER_HPP
#define NOSH_NONLINEAR_SOLVER_HPP

#include <Piro_ObserverBase.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Thyra_ModelEvaluatorDefaultBase.hpp>
#include <Tpetra_Vector.hpp>

#include <boost/any.hpp>

#include <map>
#include <memory>

using list = std::map<std::string, boost::any>;
namespace nosh {
  class nonlinear_observer: public Piro::ObserverBase<double>
  {
    public:
    nonlinear_observer(): solution(Teuchos::null) {};

    virtual ~nonlinear_observer() {};

    using Piro::ObserverBase<double>::observeSolution;
    virtual
    void
    observeSolution(const Thyra::VectorBase<double> &sol)
    {
      this->solution = sol.clone_v();
    }

    public:
    Teuchos::RCP<const Thyra::VectorBase<double>> solution;
  };

  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  nonlinear_solve(
      const std::shared_ptr<Thyra::ModelEvaluatorDefaultBase<double>> & model,
      std::map<std::string, boost::any> solver_params
      );

  void
  parameter_continuation(
      const std::shared_ptr<Thyra::ModelEvaluatorDefaultBase<double>> & model,
      std::map<std::string, boost::any> solver_params
      );
}
#endif // NOSH_NONLINEAR_SOLVER_HPP
