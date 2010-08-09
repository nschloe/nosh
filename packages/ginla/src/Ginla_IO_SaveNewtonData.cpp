#include "Ginla_IO_SaveNewtonData.h"

#include "Ginla_IO_StateWriter.h"
#include "Ginla_StateTranslator.h"

#include <NOX_Solver_Generic.H>
#include <NOX_Epetra_Group.H>

// =============================================================================
Ginla::IO::SaveNewtonData::
SaveNewtonData ( const Teuchos::RCP<const Ginla::IO::StateWriter> & stateWriter,
                 const Teuchos::RCP<const Ginla::StateTranslator> & translator
               ) :
        stateWriter_ ( stateWriter ),
        translator_ ( translator )
{
}
// =============================================================================
Ginla::IO::SaveNewtonData::
~SaveNewtonData()
{
}
// =============================================================================
void
Ginla::IO::SaveNewtonData::
runPostIterate ( const NOX::Solver::Generic& solver )
{
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& solGrp =
        dynamic_cast<const NOX::Epetra::Group&> ( solver.getSolutionGroup() );

    const Epetra_Vector& currentSol =
        ( dynamic_cast<const NOX::Epetra::Vector&> ( solGrp.getX() ) ).getEpetraVector();

    stateWriter_->write( translator_->createState(currentSol),
                         solver.getNumIterations() );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    TEST_FOR_EXCEPTION ( !solGrp.isF(),
                         std::logic_error,
                         "Group contains invalid right hand side. Has F ever been computed?" );
    const Epetra_Vector& currentResidual =
        ( dynamic_cast<const NOX::Epetra::Vector&> ( solGrp.getF() ) ).getEpetraVector();

    stateWriter_->write( translator_->createState(currentResidual),
                         solver.getNumIterations(),
                         "-residual" );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return;
}
// =============================================================================
