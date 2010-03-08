#include "GL_IO_SaveNewtonData.h"

#include <Teuchos_ParameterList.hpp>
#include <NOX_Solver_Generic.H>

#include <NOX_Epetra_Group.H>

// =============================================================================
GL::IO::SaveNewtonData::
SaveNewtonData ( const Teuchos::RCP<const AbstractStateWriter>  & stateWriter,
                 const std::string                              & outputDir
               ) :
        numRunPreIterate ( 0 ),
        stateWriter_ ( stateWriter ),
        outputDir_ ( outputDir )
{
}
// =============================================================================
GL::IO::SaveNewtonData::
~SaveNewtonData()
{
}
// =============================================================================
void GL::IO::SaveNewtonData::
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
    << outputDir_ << "/newton-sol-" << numRunPreIterate;
    stateWriter_->writeSolutionToFile ( currentSol,
                                        fileStream.str() );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    TEST_FOR_EXCEPTION ( !solGrp.isF(),
                         std::logic_error,
                         "Group contains invalid right hand side. Has F ever been computed?" );
    const Epetra_Vector& currentResidual =
        ( dynamic_cast<const NOX::Epetra::Vector&> ( solGrp.getF() ) ).getEpetraVector();
    fileStream.str(std::string()); // empty the filestream
    fileStream
    << outputDir_ << "/newton-res-" << numRunPreIterate;
    stateWriter_->writeAbstractStateToFile ( currentResidual,
                                             fileStream.str() );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
// =============================================================================
