#include "LinearSolver.hpp"

#include "LinearOperator.hpp"

#include <Teuchos_RCP.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_ParameterList.hpp>

#include <map>

// =============================================================================
void
Nosh::
linearSolve(
    const Nosh::LinearOperator & A,
    const Nosh::Expression & f,
    Nosh::Function & x,
    std::map<std::string, boost::any> solverParams
    )
{
  // create f vector
  auto fVec = Nosh::integrateOverControlVolumes(f, *A.mesh);

  // apply boundary conditions to f
  const auto boundaryNodes = A.mesh->getBoundaryNodes();
  const VectorFieldType & coordsField = A.mesh->getNodeField("coordinates");
  for (const auto boundaryNode: boundaryNodes) {
    const auto coord = A.mesh->getNodeValue(coordsField, boundaryNode);
    for (const auto & bc: A.bcs) {
      if (bc->isInside(coord)) {
        const auto gid = A.mesh->gid(boundaryNode);
        fVec->replaceGlobalValue(gid, bc->eval(coord));
        break; // only set one bc per boundary point
      }
    }
  }

  //auto out = Teuchos::VerboseObjectBase::getDefaultOStream();
  //A.describe(*out, Teuchos::VERB_EXTREME);
  //return;

  Stratimikos::DefaultLinearSolverBuilder builder;
  auto p = Teuchos::rcp(new Teuchos::ParameterList());
  stdmap2teuchoslist(solverParams, *p);
  builder.setParameterList(p);
  auto lowsFactory = builder.createLinearSolveStrategy("");
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  const Tpetra::Operator<double,int,int> & opA = A;
  auto thyraA = Thyra::createConstLinearOp(Teuchos::rcpFromRef(opA)); // throws
  auto lows = Thyra::linearOpWithSolve(
      *lowsFactory,
      thyraA
      );
  const Tpetra::Vector<double,int,int> & vecF = *fVec;
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
std::map<std::string, boost::any>
Nosh::
defaultLinearSolverParams()
{
  using list = std::map<std::string, boost::any>;
  return {
    {"Linear Solver Type", "Belos"},
    {"Linear Solver Types", list{
      {"Belos", list{
        {"Solver Type", "Pseudo Block GMRES"},
        {"Solver Types", list{
          {"Pseudo Block GMRES", list{
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
Nosh::
stdmap2teuchoslist(
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
        stdmap2teuchoslist(
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
