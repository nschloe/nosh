#ifndef NOSH_LINEARSOLVER_HPP
#define NOSH_LINEARSOLVER_HPP

#include "LinearOperator.hpp"
#include "Function.hpp"

namespace Nosh {

  Teuchos::RCP<Teuchos::ParameterList>
  defaultLinearSolverParams();

  void
  linearSolve(
      Nosh::LinearOperator & A,
      Nosh::Function & f,
      Nosh::Function & x,
      Teuchos::RCP<Teuchos::ParameterList> solverParams = Nosh::defaultLinearSolverParams()
      );
}
#endif // NOSH_LINEARSOLVER_HPP
