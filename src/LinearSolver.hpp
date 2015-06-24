#ifndef NOSH_LINEARSOLVER_HPP
#define NOSH_LINEARSOLVER_HPP

#include "LinearOperator.hpp"
#include "Expression.hpp"
#include "Function.hpp"

namespace Nosh {

  Teuchos::RCP<Teuchos::ParameterList>
  defaultLinearSolverParams();

  void
  linearSolve(
      const Nosh::LinearOperator & A,
      const Nosh::Expression & f,
      Nosh::Function & x,
      Teuchos::RCP<Teuchos::ParameterList> solverParams = Nosh::defaultLinearSolverParams()
      );
}
#endif // NOSH_LINEARSOLVER_HPP
