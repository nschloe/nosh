#ifndef NOSH_LINEARSOLVER_HPP
#define NOSH_LINEARSOLVER_HPP

#include "linear_operator.hpp"
#include "expression.hpp"
#include "function.hpp"

#include <Teuchos_ParameterList.hpp>
#include <boost/any.hpp>

#include <map>

namespace nosh {
  std::map<std::string, boost::any> default_linear_solver_params();

  void
  linear_solve(
      const nosh::linear_operator & A,
      const nosh::expression & f,
      nosh::function & x,
      std::map<std::string, boost::any> solver_params = nosh::default_linear_solver_params()
      );

  void
  std_map_to_teuchos_list(
      const std::map<std::string, boost::any> & map,
      Teuchos::ParameterList & p
      );
}
#endif // NOSH_LINEARSOLVER_HPP
