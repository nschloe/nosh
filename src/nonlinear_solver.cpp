#include "nonlinear_solver.hpp"
#include "helper.hpp"

#include <Piro_NOXSolver.hpp>
#include <Piro_LOCASolver.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <map>


// =============================================================================
void
nosh::
nonlinear_solve(
    const std::shared_ptr<Thyra::ModelEvaluatorDefaultBase<double>> & model,
    std::map<std::string, boost::any> solver_params
    )
{
  auto p = Teuchos::rcp(new Teuchos::ParameterList());
  std_map_to_teuchos_list(solver_params, *p);

  auto piro = std::make_shared<Piro::NOXSolver<double>>(
        p,
        Teuchos::rcp(model)
        // Teuchos::rcp(observer)
        );

  // Now the setting of inputs and outputs.
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = piro->createInArgs();
  inArgs.set_p(
      0,
      piro->getNominalValues().get_p(0)
      );

  // Set output arguments to evalModel call.
  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = piro->createOutArgs();

  // Now solve the problem and return the responses.
  piro->evalModel(inArgs, outArgs);

  return;
}
// =============================================================================
