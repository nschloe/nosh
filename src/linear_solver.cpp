#include "linear_solver.hpp"

#include "matrix.hpp"

#include <Teuchos_RCP.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_ParameterList.hpp>

#include <map>

// =============================================================================
void
nosh::
linear_solve(
    const nosh::matrix & A,
    const nosh::expression & f,
    nosh::function & x,
    std::map<std::string, boost::any> solver_params
    )
{
  // create f vector
  auto f_vec = nosh::integrate_over_control_volumes(f, *A.mesh);
  // solve
  linear_solve(A, f_vec, x, solver_params);
  return;
}
// =============================================================================
void
nosh::
linear_solve(
    const nosh::matrix & A,
    std::shared_ptr<Tpetra::Vector<double,int,int>> f_vec,
    nosh::function & x,
    std::map<std::string, boost::any> solver_params
    )
{
  // apply boundary conditions to f
  const auto boundary_nodes = A.mesh->boundary_nodes();
  const vector_fieldType & coords_field = A.mesh->get_node_field("coordinates");
  for (const auto boundary_node: boundary_nodes) {
    const auto coord = A.mesh->get_node_value(coords_field, boundary_node);
    for (const auto & bc: A.bcs) {
      if (bc->is_inside(coord)) {
        const auto gid = A.mesh->gid(boundary_node);
        f_vec->replaceGlobalValue(gid, bc->eval(coord));
        break; // only set one bc per boundary point
      }
    }
  }

  //auto out = Teuchos::VerboseObjectBase::getDefaultOStream();
  //A.describe(*out, Teuchos::VERB_EXTREME);
  //return;

  // set x to 0
  x.putScalar(0.0);

  Stratimikos::DefaultLinearSolverBuilder builder;
  auto p = Teuchos::rcp(new Teuchos::ParameterList());
  std_map_to_teuchos_list(solver_params, *p);
  builder.setParameterList(p);
  auto lowsFactory = builder.createLinearSolveStrategy("");
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  const Tpetra::Operator<double,int,int> & opA = A;
  auto thyraA = Thyra::createConstLinearOp(Teuchos::rcpFromRef(opA)); // throws
  auto lows = Thyra::linearOpWithSolve(
      *lowsFactory,
      thyraA
      );
  const Tpetra::Vector<double,int,int> & vecF = *f_vec;
  Tpetra::Vector<double,int,int> & vecX = x;
  auto status = Thyra::solve<double>(
      *lows,
      Thyra::NOTRANS,
      *Thyra::createConstVector(Teuchos::rcpFromRef(vecF)),
      Thyra::createVector(Teuchos::rcpFromRef(vecX)).ptr()
      );
  std::cout << "Solve status: " << status << std::endl;
  return;
}
// =============================================================================
void
nosh::
scaled_linear_solve(
    nosh::matrix & A,
    const nosh::expression & f,
    nosh::function & x,
    std::map<std::string, boost::any> solver_params
    )
{
  // create f vector
  auto b = nosh::integrate_over_control_volumes(f, *A.mesh);

  // create scaling vectors
  const auto control_volumes = A.mesh->control_volumes();
  auto inv_sqrt_control_volumes =
    Tpetra::Vector<double,int,int>(control_volumes->getMap());

  auto cv_data = control_volumes->getData();
  auto inv_sqrt_cv_data = inv_sqrt_control_volumes.getDataNonConst();
  for (int k = 0; k < cv_data.size(); k++) {
    inv_sqrt_cv_data[k] = 1.0 / sqrt(cv_data[k]);
  }

  // scale A
  A.leftScale(inv_sqrt_control_volumes);
  A.rightScale(inv_sqrt_control_volumes);

  // scale fvec
  auto b_data = b->getDataNonConst();
  for (int k = 0; k < b_data.size(); k++) {
    b_data[k] *= inv_sqrt_cv_data[k];
  }

  // solve
  linear_solve(A, b, x, solver_params);

  // scale the solution
  auto x_data = x.getDataNonConst();
  for (int k = 0; k < x_data.size(); k++) {
    x_data[k] *= inv_sqrt_cv_data[k];
  }

  return;
}
// =============================================================================
std::map<std::string, boost::any>
nosh::
default_linear_solver_params()
{
  using list = std::map<std::string, boost::any>;
  return {
    {"Linear Solver Type", "Belos"},
    {"Linear Solver Types", list{
      {"Belos", list{
        {"Solver Type", "Pseudo Block GMRES"},
        {"Solver Types", list{
          {"Pseudo Block GMRES", list{
            {"Convergence Tolerance", 1.0e-10},
            {"Output Frequency", 1},
            {"Output Style", 1},
            {"Verbosity", 33}
          }}
        }}
      }}
    }},
    {"Preconditioner Type", "None"}
  };

  // auto p = Teuchos::rcp(new Teuchos::ParameterList);
  // p->set("Linear Solver Type", "Belos");
  // auto & belosList =
  //   p->sublist("Linear Solver Types")
  //   .sublist("Belos");
  // belosList.set("Solver Type", "Pseudo Block GMRES");

  // auto & solverList =
  //   belosList.sublist("Solver Types")
  //   .sublist("Pseudo Block CG");
  // solverList.set("Output Frequency", 1);
  // solverList.set("Output Style", 1);
  // solverList.set("Verbosity", 33);

  // p->set("Preconditioner Type", "None");

  // return p;
}
// =============================================================================
void
nosh::
std_map_to_teuchos_list(
    const std::map<std::string, boost::any> & map,
    Teuchos::ParameterList & p
    )
{
  for (const auto & entry: map) {
      if(entry.second.type() == typeid(int)) {
        p.set(entry.first, boost::any_cast<int>(entry.second));
      } else if(entry.second.type() == typeid(double)) {
        p.set(entry.first, boost::any_cast<double>(entry.second));
      } else if(entry.second.type() == typeid(const char*)) {
        p.set(entry.first, boost::any_cast<const char*>(entry.second));
      } else if(entry.second.type() == typeid(std::string)) {
        p.set(entry.first, boost::any_cast<std::string>(entry.second));
      } else if(entry.second.type() == typeid(std::map<std::string, boost::any>)) {
        std_map_to_teuchos_list(
            boost::any_cast<std::map<std::string, boost::any>>(entry.second),
            p.sublist(entry.first)
            );
      } else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            true,
            "Unknown value type of key \"" << entry.first << "\"."
            );
      }
  }
  return;
}
// =============================================================================
