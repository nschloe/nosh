#ifndef NOSH_LINEARSOLVER_HPP
#define NOSH_LINEARSOLVER_HPP

#include "LinearOperator.hpp"
#include "Expression.hpp"
#include "Function.hpp"

#include <Teuchos_ParameterList.hpp>
#include <boost/any.hpp>

#include <map>

namespace Nosh {
  std::map<std::string, boost::any>
  defaultLinearSolverParams();

  void
  linearSolve(
      const Nosh::LinearOperator & A,
      const Nosh::Expression & f,
      Nosh::Function & x,
      std::map<std::string, boost::any> solverParams = Nosh::defaultLinearSolverParams()
      );

  void
  stdmap2teuchoslist(
      const std::map<std::string, boost::any> & map,
      Teuchos::ParameterList & p
      );
}
#endif // NOSH_LINEARSOLVER_HPP
