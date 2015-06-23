#include "LinearSolver.hpp"

#include "LinearOperator.hpp"

#include <Teuchos_RCP.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Tpetra_CrsMatrix.hpp>

// =============================================================================
void
Nosh::
linearSolve(
    const Nosh::LinearOperator & A,
    Nosh::Function & f,
    Nosh::Function & x,
    Teuchos::RCP<Teuchos::ParameterList> solverParams
    )
{
  // apply boundary conditions to f
  const auto boundaryNodes = A.mesh->getBoundaryNodes();
  const VectorFieldType & coordsField = A.mesh->getNodeField("coordinates");
  for (const auto boundaryNode: boundaryNodes) {
    const auto coord = A.mesh->getNodeValue(coordsField, boundaryNode);
    for (const auto & bc: A.bcs) {
      if (bc->isInside(coord)) {
        const auto gid = A.mesh->gid(boundaryNode);
        f.replaceGlobalValue(gid, bc->eval(coord));
        break; // only set one bc per boundary point
      }
    }
  }

  //auto out = Teuchos::VerboseObjectBase::getDefaultOStream();
  //A.describe(*out, Teuchos::VERB_EXTREME);
  //return;

  Stratimikos::DefaultLinearSolverBuilder builder;
  builder.setParameterList(solverParams);
  auto lowsFactory = builder.createLinearSolveStrategy("");
  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  const Tpetra::Operator<double,int,int> & opA = A;
  auto thyraA = Thyra::createConstLinearOp(Teuchos::rcpFromRef(opA)); // throws
  auto lows = Thyra::linearOpWithSolve(
      *lowsFactory,
      thyraA
      );
  const Tpetra::Vector<double,int,int> & vecF = f;
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
Teuchos::RCP<Teuchos::ParameterList>
Nosh::
defaultLinearSolverParams()
{
  auto p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "Belos");
  auto & belosList =
    p->sublist("Linear Solver Types")
    .sublist("Belos");
  belosList.set("Solver Type", "Pseudo Block GMRES");

  auto & solverList =
    belosList.sublist("Solver Types")
    .sublist("Pseudo Block CG");
  solverList.set("Output Frequency", 1);
  solverList.set("Output Style", 1);
  solverList.set("Verbosity", 33);

  p->set("Preconditioner Type", "None");

  return p;
}
// =============================================================================
