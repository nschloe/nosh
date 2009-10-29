#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>

#include <EpetraExt_RowMatrixOut.h>

#include "glNox.h"

int main ( int argc, char *argv[] )
{

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init ( &argc,&argv );
#endif

  // Create a communicator for Tpetra objects
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm
                                         = Teuchos::DefaultComm<int>::getComm();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_MpiComm> eComm
  = Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
  Teuchos::RCP<Epetra_SerialComm>  eComm
  = Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

  // ===========================================================================
  // handle command line arguments
  Teuchos::CommandLineProcessor My_CLP;

  My_CLP.setDocString (
    "This program solves the Ginzburg--Landau problem with a NOX interace.\n"
  );

  bool verbose = false;
  My_CLP.setOption ( "verbose", "silent",
                     &verbose,
                     "Verbostity flag" );

  bool matlabMatrix = false;
  My_CLP.setOption ( "jacobian-file", "no-jacobian-file",
                     &matlabMatrix,
                     "Save the jacobian in a text file" );

  bool computeEigenvalues = false;
  My_CLP.setOption ( "eigenvalues", "no-eigenvalues",
                     &computeEigenvalues,
                     "Compute eigenvalue approximations in the solution" );
  
  bool computeConditionNumber = false;
  My_CLP.setOption ( "condest", "no-condest",
                     &computeConditionNumber,
                     "Compute condition number approximations in the solution" );

  std::string filename = "";
  My_CLP.setOption ( "input-guess",
                     &filename,
                     "File name with initial guess" );

  std::string outputdir = "data";
  My_CLP.setOption ( "output-dir",
                     &outputdir,
                     "Directory to which the solution files are written" );

  // print warning for unrecognized arguments
  My_CLP.recogniseAllOptions ( true );

  // don't throw exceptions
  My_CLP.throwExceptions ( false );

  // finally, parse the stuff!
  Teuchos::CommandLineProcessor::EParseCommandLineReturn
  parseReturn = My_CLP.parse ( argc, argv );
  if ( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
  }

  if ( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 1; // Error!;
  }
  // ===========================================================================


  // create GL-NOX object with initial parameters/guess
  Teuchos::RCP<glNox> myNoxObject = Teuchos::ENull();
  if ( filename.length()>0 ) {
      myNoxObject = Teuchos::rcp( new glNox( filename, Comm, eComm ) );
  }
  else {
      int    Nx         = 50;
      double edgeLength = 10.0;
      double H0         = 0.4;
      myNoxObject = Teuchos::rcp( new glNox( Nx, edgeLength, H0, Comm, eComm ) );
  }

  myNoxObject->setVerbose( verbose );

  // set default solver options
  int maxNonlinearIterations = 30;
  myNoxObject->setSolverOptions( maxNonlinearIterations );

  myNoxObject->createSolverGroup();
  myNoxObject->createConvergenceTests();
  myNoxObject->createSolver();

  // solve the system
  myNoxObject->solve();

  // compute the condition number
  if ( computeConditionNumber ) {
      double kappa = myNoxObject->computeJacobianConditionNumber();
      std::cout << "Condition number: kappa = " << kappa << "." << std::endl;
  }

  // spit out the Jacobian in MATLAB readable format
  if ( matlabMatrix ) {
      std::string jacFilename = "jacobianMatrix.dat";
//       EpetraExt::RowMatrixToMatlabFile(jacFilename.c_str(),*(glsystem->getJacobian()));
  }

  // compute the eigenvaluesof the Jacobian
  if ( computeEigenvalues )
        myNoxObject->computeJacobianEigenvalues();

  
  // print the solution to a file
  myNoxObject->printSolutionToFile( "data/solution.vtk" );

  // check the convergence status
  int status = myNoxObject->checkConvergence();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return status;
}