#include "linear_solver.hpp"

#include "matrix.hpp"

#include <Amesos2.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>

//#include <Amesos2_Details_LinearSolverFactory.hpp>
//#include <Belos_Details_LinearSolverFactory.hpp>
//#include <Ifpack2_Details_LinearSolverFactory.hpp>
// #include <Trilinos_Details_LinearSolver.hpp>
// #include <Trilinos_Details_LinearSolverFactory.hpp>

#include <map>

using SC = double;
using LO = int;
using GO = int;
using MV = Tpetra::MultiVector<SC,LO,GO>;
using OP = Tpetra::CrsMatrix<SC,LO,GO>;
using NormType = MV::mag_type;


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
  auto b = nosh::integrate_over_control_volumes(f, *A.mesh);
  // solve
  linear_solve(A, b, x, solver_params);
  return;
}
// =============================================================================
void
nosh::
linear_solve(
    const nosh::matrix & A,
    std::shared_ptr<Tpetra::Vector<double,int,int>> b,
    nosh::function & x,
    std::map<std::string, boost::any> solver_params
    )
{
  // apply boundary conditions to b
  const auto boundary_nodes = A.mesh->boundary_nodes();
  for (const auto boundary_node: boundary_nodes) {
    const auto coord = A.mesh->get_coords(boundary_node);
    for (const auto & bc: A.bcs) {
      TEUCHOS_ASSERT(bc != nullptr);
      if (bc->is_inside(coord)) {
        const auto gid = A.mesh->gid(boundary_node);
        // TODO don't check here but only get the array of owned boundary nodes
        // in the first place
        if (b->getMap()->isNodeGlobalElement(gid)) {
          b->replaceGlobalValue(gid, bc->eval(coord));
          break; // only set one bc per boundary point
        }
      }
    }
  }

  //auto A_rcp = Teuchos::rcpFromRef(A);
  //solver->setMatrix(A_rcp);
  //// solver.setParameters(); // TODO
  //solver->solve(x, *b);

  const std::string package =
    boost::any_cast<const char *>(solver_params.at("package"));
  if (package == "Amesos2") {
      linear_solve_amesos2(A, b, x, solver_params);
  } else if (package == "Belos") {
      linear_solve_belos(A, b, x, solver_params);
  } else if (package == "MueLu") {
      linear_solve_muelu(A, b, x, solver_params);
  } else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          "Unknown linear solver package \"" << package << "\"."
          );
  }
  return;
}
// =============================================================================
void
nosh::
linear_solve_amesos2(
    const nosh::matrix & A,
    std::shared_ptr<Tpetra::Vector<double,int,int>> b,
    nosh::function & x,
    std::map<std::string, boost::any> solver_params
    )
{
  auto solver = Amesos2::create<OP,MV>(
        "Superlu",
        Teuchos::rcpFromRef(A),
        Teuchos::rcpFromRef(x),
        Teuchos::rcp(b)
        );

  //// Create a Teuchos::ParameterList to hold solver parameters
  //Teuchos::ParameterList amesos2_params("Amesos2");
  //Teuchos::ParameterList superlu_params = amesos2_params.sublist("SuperLU");
  //superlu_params.set("Trans","TRANS","Whether to solve with A^T");
  //superlu_params.set("Equil",false,"Whether to equilibrate the system before solve");
  //superlu_params.set("ColPerm","NATURAL","Use 'natural' ordering of columns");
  //solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );

  solver->symbolicFactorization().numericFactorization().solve();

  return;
}
// =============================================================================
void
nosh::
linear_solve_belos(
    const nosh::matrix & A,
    std::shared_ptr<Tpetra::Vector<double,int,int>> b,
    nosh::function & x,
    std::map<std::string, boost::any> solver_params
    )
{
  // set x to 0
  x.putScalar(0.0);

  Stratimikos::DefaultLinearSolverBuilder builder;
  auto p = Teuchos::rcp(new Teuchos::ParameterList());
  std_map_to_teuchos_list(convert_to_belos_parameters(solver_params), *p);
  builder.setParameterList(p);
  auto lowsFactory = builder.createLinearSolveStrategy("");
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  const Tpetra::Operator<double,int,int> & opA = A;
  auto thyraA = Thyra::createConstLinearOp(Teuchos::rcpFromRef(opA)); // throws
  auto lows = Thyra::linearOpWithSolve(
      *lowsFactory,
      thyraA
      );
  const Tpetra::Vector<double,int,int> & vecF = *b;
  Tpetra::Vector<double,int,int> & vecX = x;
  auto status = Thyra::solve<double>(
      *lows,
      Thyra::NOTRANS,
      *Thyra::createConstVector(Teuchos::rcpFromRef(vecF)),
      Thyra::createVector(Teuchos::rcpFromRef(vecX)).ptr()
      );

  if (A.getComm()->getRank() == 0) {
    std::cout << status << std::endl;
  }
  return;
}
// =============================================================================
std::shared_ptr<MueLu::Hierarchy<double,int,int>>
nosh::
get_muelu_hierarchy(
    const nosh::matrix & A,
    const std::map<std::string, boost::any> & muelu_params
    )
{
  // Tpetra -> Xpetra
  Teuchos::RCP<const Tpetra::CrsMatrix<double,int,int>> ATpetra =
    Teuchos::rcpFromRef(A);
  // cast away the const from A :(
  auto nonconst_ATpetra =
    Teuchos::rcp_const_cast<Tpetra::CrsMatrix<double,int,int>>(ATpetra);
  auto AXpetra = MueLu::TpetraCrs_To_XpetraMatrix(nonconst_ATpetra);

  auto map = AXpetra->getRowMap();

  auto p = Teuchos::rcp(new Teuchos::ParameterList());
  const auto & params = boost::any_cast<std::map<std::string,boost::any>>(
      muelu_params
      );
  std_map_to_teuchos_list(params, *p);
  std::cout << p << std::endl;

  auto mueLuFactory =
    MueLu::ParameterListInterpreter<double,int,int>(*p, map->getComm());

  auto H = Teuchos::get_shared_ptr(mueLuFactory.CreateHierarchy());
  H->GetLevel(0)->Set("A", AXpetra);

  //// build null space vector
  //auto nullspace = Xpetra::MultiVectorFactory<double,int,int>::Build(map, 1);
  //nullspace->putScalar(1.0);
  //H->GetLevel(0)->Set("Nullspace", nullspace);

  //// TODO
  //// get the coordinates as multivector
  //RCP<MultiVector> coords = Teuchos::rcp(new Xpetra::EpetraMultiVector(epCoord));
  //H->GetLevel(0)->Set("Coordinates", coords);

  mueLuFactory.SetupHierarchy(*H);

  return H;
}
// =============================================================================
void
nosh::
linear_solve_muelu(
    const nosh::matrix & A,
    std::shared_ptr<Tpetra::Vector<double,int,int>> b,
    nosh::function & x,
    std::map<std::string, boost::any> solver_params
    )
{
  x.putScalar(0.0);
  // Tpetra -> Xpetra
  auto bXpetra = Xpetra::toXpetra(Teuchos::rcpFromRef(*b));
  Tpetra::Vector<double,int,int> & xTpetra = x;
  auto xXpetra = Xpetra::toXpetra(Teuchos::rcpFromRef(xTpetra));

  std::map<std::string, boost::any> muelu_params;
  try {
    muelu_params = boost::any_cast<std::map<std::string, boost::any>>(
        solver_params.at("parameters")
        );
  } catch (std::out_of_range) {
    muelu_params = {};
  }
  auto H = get_muelu_hierarchy(A, muelu_params);
  H->IsPreconditioner(false);

  const int mgridSweeps = 100;
  H->Iterate(*bXpetra, *xXpetra, mgridSweeps);

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
convert_to_belos_parameters(
    const std::map<std::string, boost::any> & in_map
    )
{
  if (in_map.find("method") != in_map.end()) {
    const std::string method =
      boost::any_cast<const char *>(in_map.at("method"));
    std::map<std::string, boost::any> out_map = {
      {"Linear Solver Type", "Belos"},
      {"Linear Solver Types", list{
        {"Belos", list{
          {"Solver Type", method}
        }}
      }},
      {"Preconditioner Type", "None"}
    };

    auto lst = boost::any_cast<list>(out_map.at("Linear Solver Types"));
    auto belos = boost::any_cast<list>(lst.at("Belos"));
    if (in_map.find("parameters") != in_map.end()) {
      belos.insert({
          "Solver Types",
          list{{method, in_map.at("parameters")}}
          });
    } else {
      // insert default parameters
      belos.insert({
          "Solver Types",
          list{{method,
            list{
              {"Convergence Tolerance", 1.0e-10},
              {"Output Frequency", 1},
              {"Output Style", 1},
              {"Verbosity", 33}
            }
            }}
          });
    }
    return out_map;
  } else {
    return {};
  }
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
