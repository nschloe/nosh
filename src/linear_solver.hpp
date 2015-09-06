#ifndef NOSH_LINEARSOLVER_HPP
#define NOSH_LINEARSOLVER_HPP

#include "matrix.hpp"
#include "expression.hpp"
#include "function.hpp"

#include <Teuchos_ParameterList.hpp>
#include <boost/any.hpp>

#include <map>

using list = std::map<std::string, boost::any>;
namespace nosh {
  // https://trilinos.org/docs/dev/packages/stratimikos/doc/html/index.html
  // http://stackoverflow.com/a/14425299/353337
  static
  const std::map<std::string, boost::any> default_linear_solver_params = {
    {"method", "Pseudo Block GMRES"},
    {"parameters", list{
      {"Convergence Tolerance", 1.0e-10},
      {"Output Frequency", 1},
      {"Output Style", 1},
      {"Verbosity", 33}
      }}
  };
  //std::map<std::string, boost::any> default_linear_solver_params =
  //{
  //  {"Linear Solver Type", "Belos"},
  //  {"Linear Solver Types", list{
  //    {"Belos", list{
  //      {"Solver Type", "Pseudo Block GMRES"},
  //      {"Solver Types", list{
  //        {"Pseudo Block GMRES", list{
  //          {"Convergence Tolerance", 1.0e-10},
  //          {"Output Frequency", 1},
  //          {"Output Style", 1},
  //          {"Verbosity", 33}
  //        }}
  //      }}
  //    }}
  //  }},
  //  {"Preconditioner Type", "None"}
  //};

  void
  linear_solve(
      const nosh::matrix & A,
      const nosh::expression & f,
      nosh::function & x,
      std::map<std::string, boost::any> solver_params = nosh::default_linear_solver_params
      );

  void
  linear_solve(
      const nosh::matrix & A,
      std::shared_ptr<Tpetra::Vector<double,int,int>> f_vec,
      nosh::function & x,
      std::map<std::string, boost::any> solver_params = nosh::default_linear_solver_params
      );

  void
  scaled_linear_solve(
      nosh::matrix & A,
      const nosh::expression & f,
      nosh::function & x,
      std::map<std::string, boost::any> solver_params = nosh::default_linear_solver_params
      );

  std::map<std::string, boost::any>
  convert_to_belos_parameters(
      const std::map<std::string, boost::any> & map
      );

  void
  std_map_to_teuchos_list(
      const std::map<std::string, boost::any> & map,
      Teuchos::ParameterList & p
      );
}
#endif // NOSH_LINEARSOLVER_HPP
