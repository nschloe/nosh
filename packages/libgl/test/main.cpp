#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <EpetraExt_RowMatrixOut.h>

// =============================================================================
int main ( int argc, char *argv[] )
{ 
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( &argc,&argv );
#endif

    // Create a communicator for Tpetra objects
    const Teuchos::RCP<const Teuchos::Comm<int> > Comm =
    Teuchos::DefaultComm<int>::getComm();

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm =
    Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm>  eComm =
    Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

    int status = 0;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return status;
}
// =========================================================================
