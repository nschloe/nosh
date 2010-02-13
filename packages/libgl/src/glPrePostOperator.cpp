#include "glPrePostOperator.h"

#include <Teuchos_ParameterList.hpp>
#include <NOX_Solver_Generic.H>

#include <NOX_Epetra_Group.H>

// =============================================================================
GlPrePostOperator::GlPrePostOperator ( const Teuchos::RCP<const AbstractStateWriter>  & stateWriter,
                                       const std::string                              & outputDir,
                                       const std::string                              & outputFormat ) :
        numRunPreIterate ( 0 ),
        stateWriter_ ( stateWriter ),
        outputDir_ ( outputDir ),
        filenameExtension_ ( "" )
{
  if ( outputFormat.compare("VTK")==0 )
    filenameExtension_ = "vtk";
  else if ( outputFormat.compare("VTI")==0 )
    filenameExtension_ = "vti";
  else
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Illegal output format \"" << outputFormat << "\"." );
}
// =============================================================================
GlPrePostOperator::~GlPrePostOperator()
{
}
// =============================================================================
void GlPrePostOperator::
runPostIterate ( const NOX::Solver::Generic& solver )
{
    string fileName;

    ++numRunPreIterate;

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& solGrp =
        dynamic_cast<const NOX::Epetra::Group&> ( solver.getSolutionGroup() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    const Epetra_Vector& currentSol =
        ( dynamic_cast<const NOX::Epetra::Vector&> ( solGrp.getX() ) ).getEpetraVector();
    stringstream fileStream;
    fileStream
    << outputDir_ << "/newton-sol-" << numRunPreIterate
    << "." << filenameExtension_;
    stateWriter_->writeSolutionToFile ( currentSol,
                                        fileStream.str() );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    TEST_FOR_EXCEPTION ( !solGrp.isF(),
                         std::logic_error,
                         "Group contains invalid right hand side. Has F ever been computed?" );
    const Epetra_Vector& currentResidual =
        ( dynamic_cast<const NOX::Epetra::Vector&> ( solGrp.getF() ) ).getEpetraVector();
    stringstream fileStream;
    fileStream
    << outputDir_ << "/newton-res-" << numRunPreIterate 
    << "." << filenameExtension_;
    stateWriter_->writeAbstractStateToFile ( currentResidual,
                                             fileStream.str() );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
// =============================================================================
