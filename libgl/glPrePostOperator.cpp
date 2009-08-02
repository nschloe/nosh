#include "glPrePostOperator.h"

#include "glSystem.h"

#include <Teuchos_ParameterList.hpp>
#include <NOX_Solver_Generic.H>

#include <NOX_Epetra_Group.H>

#include <EpetraExt_Utils.h> // for toString

GlPrePostOperator::GlPrePostOperator( Teuchos::RCP<GlSystem>        glsystem,
                                      const Teuchos::ParameterList& problemParams) :
  numRunPreIterate(0),
  problemParameters_( problemParams ),
  glsystem_(glsystem)
{
}

GlPrePostOperator::~GlPrePostOperator()
{
}

void GlPrePostOperator::
runPreIterate(const NOX::Solver::Generic& solver)
{
  ++numRunPreIterate;

  string fileName = "data/newton-step-"+EpetraExt::toString(numRunPreIterate)+".vtk";

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& solGrp =
    dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& currentSol =
    (dynamic_cast<const NOX::Epetra::Vector&>(solGrp.getX())).
    getEpetraVector();

  glsystem_->solutionToFile( currentSol,
                             problemParameters_,
                             fileName );

}