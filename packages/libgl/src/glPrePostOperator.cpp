#include "glPrePostOperator.h"

#include <Teuchos_ParameterList.hpp>
#include <NOX_Solver_Generic.H>

#include <NOX_Epetra_Group.H>

#include <EpetraExt_Utils.h> // for toString

// =============================================================================
GlPrePostOperator::GlPrePostOperator( const Teuchos::RCP<const AbstractStateWriter>  & stateWriter,
                                      const std::string                              & outputDir) :
  numRunPreIterate(0),
  stateWriter_(stateWriter),
  outputDir_(outputDir)
{
}
// =============================================================================
GlPrePostOperator::~GlPrePostOperator()
{
}
// =============================================================================
void GlPrePostOperator::
runPostIterate(const NOX::Solver::Generic& solver)
{
  string fileName;

  ++numRunPreIterate;

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& solGrp =
             dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  const Epetra_Vector& currentSol =
    (dynamic_cast<const NOX::Epetra::Vector&>(solGrp.getX())).
    getEpetraVector();
  fileName = outputDir_ + "/newton-sol-"+EpetraExt::toString(numRunPreIterate)+".vti";
  stateWriter_->writeSolutionToFile( currentSol,
                                  fileName );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  TEST_FOR_EXCEPTION( !solGrp.isF(),
                      std::logic_error,
                      "Group contains invalid right hand side. Has F ever been computed?" );
  const Epetra_Vector& currentResidual =
    (dynamic_cast<const NOX::Epetra::Vector&>(solGrp.getF())).getEpetraVector();
  fileName = outputDir_ + "/newton-res-"+EpetraExt::toString(numRunPreIterate)+".vti";
  stateWriter_->writeAbstractStateToFile( currentResidual,
                                       fileName );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
// =============================================================================
