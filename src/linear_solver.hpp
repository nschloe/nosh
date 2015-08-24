#ifndef NOSH_LINEARSOLVER_HPP
#define NOSH_LINEARSOLVER_HPP

#include "matrix.hpp"
#include "expression.hpp"
#include "function.hpp"

#include <Teuchos_ParameterList.hpp>
#include <boost/any.hpp>

#include <map>

namespace nosh {
  std::map<std::string, boost::any> default_linear_solver_params();

  void
  linear_solve(
      const nosh::matrix & A,
      const nosh::expression & f,
      nosh::function & x,
      std::map<std::string, boost::any> solver_params = nosh::default_linear_solver_params()
      );

  void
  linear_solve(
      const nosh::matrix & A,
      std::shared_ptr<Tpetra::Vector<double,int,int>> f_vec,
      nosh::function & x,
      std::map<std::string, boost::any> solver_params = nosh::default_linear_solver_params()
      );

  void
  scaled_linear_solve(
      nosh::matrix & A,
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
