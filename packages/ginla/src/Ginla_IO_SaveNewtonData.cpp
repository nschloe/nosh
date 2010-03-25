#include "Ginla_IO_SaveNewtonData.h"

#include <NOX_Solver_Generic.H>
#include <NOX_Epetra_Group.H>

#include "Ginla_IO_StateWriter.h"
#include "Ginla_Komplex.h"

// =============================================================================
Ginla::IO::SaveNewtonData::
SaveNewtonData ( const Teuchos::RCP<const Ginla::IO::StateWriter> & stateWriter,
                 const Teuchos::RCP<const Recti::Grid::General>   & grid,
                 const Teuchos::RCP<const Ginla::Komplex>         & komplex
               ) :
        stateWriter_ ( stateWriter ),
        grid_ ( grid ),
        komplex_( komplex )
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

    stateWriter_->write( komplex_->real2complex(currentSol),
                         grid_,
                         solver.getNumIterations() );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    TEST_FOR_EXCEPTION ( !solGrp.isF(),
                         std::logic_error,
                         "Group contains invalid right hand side. Has F ever been computed?" );
    const Epetra_Vector& currentResidual =
        ( dynamic_cast<const NOX::Epetra::Vector&> ( solGrp.getF() ) ).getEpetraVector();

    stateWriter_->write( komplex_->real2complex(currentResidual),
                         grid_,
                         solver.getNumIterations(),
                         "-residual" );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return;
}
// =============================================================================
