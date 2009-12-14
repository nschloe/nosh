#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

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
    "This program solves the Ginzburg--Landau problem with a NOX interface.\n"
  );

//  bool verbose = false;
//  My_CLP.setOption ( "verbose", "silent",
//                     &verbose,
//                     "Verbostity flag" );

  std::string xmlInputFileName = "";
  My_CLP.setOption("xml-input-file", &xmlInputFileName,
      "XML file containing the parameter list", true);

  // print warning for unrecognized arguments
  My_CLP.recogniseAllOptions(true);

  // don't throw exceptions
  My_CLP.throwExceptions(false);

  // finally, parse the stuff!
  Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
  try {
      parseReturn = My_CLP.parse(argc, argv);
  }
  catch (std::exception &e)
  {
      std::cerr << e.what() << std::endl;
      return 1;
  }
  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
      return 0;
    }
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    {
      return 1; // Error!
    }

  Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(
      new Teuchos::ParameterList);
  std::cout << "Reading parameter list from \"" << xmlInputFileName << "\"."
      << std::endl;
  Teuchos::updateParametersFromXmlFile(xmlInputFileName, paramList.get());
  // =========================================================================
  // extract data of the parameter list
  Teuchos::ParameterList& ioList = paramList->sublist("IO", true);
  std::string inputGuessFile = ioList.get<string> ("Input guess");
  bool withInitialGuess = inputGuessFile.length() > 0;
  std::string outputDirectory = ioList.get<string> ("Output directory");

  bool computeEigenvalues =paramList->sublist("Eigenvalues",true)
                                     .get("Compute Eigenvalues", false);

  bool computeConditionNumbers = paramList->sublist("Condition Numbers",true)
                                     .get("Compute Condition Numbers", false);

  bool plotEachNewtonStep = paramList->sublist("IO",true)
                                      .get("Plot each Newton step",false);

  std::string jacFilename = paramList->sublist("IO",true)
                                      .get("Jacobian MATLAB matrix file name","");
  // =========================================================================

  // create GL-NOX object with initial parameters/guess
  Teuchos::RCP<glNox> myNoxObject = Teuchos::null;
  if ( inputGuessFile.length()>0 ) {
      myNoxObject = Teuchos::rcp( new glNox( inputGuessFile, Comm, eComm ) );
  }
  else {
      int    Nx         = 50;
      double edgeLength = 10.0;
      double H0         = 0.4;
      myNoxObject = Teuchos::rcp( new glNox( Nx, edgeLength, H0, Comm, eComm ) );
  }

  // set default solver options
  myNoxObject->setSolverOptions( plotEachNewtonStep,
                                 paramList->sublist("NOX",true),
                                 outputDirectory );

  myNoxObject->createSolverGroup();
  myNoxObject->createConvergenceTests( paramList->sublist("NOX Status Test",true) );
  myNoxObject->createSolver();

  // solve the system
  myNoxObject->solve();

  // compute the condition number
  if ( computeConditionNumbers ) {
      double kappa = myNoxObject->computeJacobianConditionNumber();
      std::cout << "Condition number: kappa = " << kappa << "." << std::endl;
  }

  // spit out the Jacobian in MATLAB readable format
  if ( !jacFilename.empty() ) {
//       EpetraExt::RowMatrixToMatlabFile(jacFilename.c_str(),*(glsystem->getJacobian()));
  }

  // compute the eigenvaluesof the Jacobian
  if ( computeEigenvalues )
        myNoxObject->computeJacobianEigenvalues();
  
  // print the solution to a file
  myNoxObject->printSolutionToFile( outputDirectory + "/solution.vtk" );

  // check the convergence status
  int status = myNoxObject->checkConvergence();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return status;
}
