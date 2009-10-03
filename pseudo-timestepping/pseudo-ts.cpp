#include <NOX.H>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <EpetraExt_Utils.h> // for toString

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

// User's application specific files
#include "glSystem.h"
#include "ioFactory.h"

#include <string>

int main(int argc, char *argv[])
{

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif


  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();

  // ===========================================================================
  // handle command line arguments
  Teuchos::CommandLineProcessor My_CLP;

  My_CLP.setDocString(
    "This program solves the Ginzburg--Landau problem with a NOX interace.\n"
  );

  bool verbose=false;
  My_CLP.setOption("verbose", "silent", &verbose, "Verbostity flag" );

  std::string filename = "";
  My_CLP.setOption("input-guess", &filename, "File name with initial guess");

  std::string outputdir = "data";
  My_CLP.setOption("output-dir", &outputdir, "Directory to which the solution files are written");

  // print warning for unrecognized arguments
  My_CLP.recogniseAllOptions(true);

  // don't throw exceptions
  My_CLP.throwExceptions(false);

  // finally, parse the stuff!
  Teuchos::CommandLineProcessor::EParseCommandLineReturn
                                       parseReturn = My_CLP.parse( argc, argv );
  if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED )
    return 0;

  if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   )
    return 1; // Error!

  bool withInitialGuess = filename.length()>0;
  // ===========================================================================

  // The following is actually only necessary when no input file is given, but
  // the VTK format as adapted here has the shortcoming that for the parameters,
  // it does not contain the data type (double, int,...). Hence, the list must
  // be present beforehand to check back for existing parameter names and their
  // types.
  Teuchos::ParameterList      problemParameters;
  problemParameters.set("Nx",50);
  problemParameters.set("edgelength",10.0);
  problemParameters.set("H0",0.4);


  // ---------------------------------------------------------------------------
  std::vector<double_complex> psiLexicographic;
  if (withInitialGuess) {
      IoVirtual* fileIo = IoFactory::createFileIo( filename );
      fileIo->read( &psiLexicographic,
                    &problemParameters );

      delete fileIo;
  }
  // ---------------------------------------------------------------------------


  // create the gl problem
  GinzburgLandau glProblem = GinzburgLandau( problemParameters.get<int>("Nx"),
                                             problemParameters.get<double>("edgelength"),
                                             problemParameters.get<double>("H0")
                                           );

  int NumUnknowns = glProblem.getStaggeredGrid()->getNumComplexUnknowns();

  // ---------------------------------------------------------------------------
  std::vector<double_complex> psi(NumUnknowns);
  if (withInitialGuess) {
      // If there was is an initial guess, make sure to get the ordering correct.
      std::vector<int> p(NumUnknowns);
      // fill p:
      glProblem.getStaggeredGrid()->lexicographic2grid( &p );
      for (int k=0; k<NumUnknowns; k++)
          psi[p[k]] = psiLexicographic[k];
  }
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // set-up the random numbers
  srand ( time(NULL) );
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // run the loop
  int    maxSteps  = 5000;
  double delta     = 0.001;
  double tol       = 1e-10;
  bool   converged = false;
  std::vector<double_complex> update(NumUnknowns);
  for( int k=0; k<maxSteps; k++ ) {

      if (verbose) {
          // print the solution to a file
          std::string fileName = outputdir+"/"+"pseudotimestep-"+EpetraExt::toString(k)+".vtk";
          IoVirtual* fileIo = IoFactory::createFileIo( fileName );
          fileIo->write( psi,
                         problemParameters,
                         *(glProblem.getStaggeredGrid()) );
      }

      // get the update
      for ( int l=0; l<NumUnknowns; l++ )
          update[l] = glProblem.computeGl( l, psi );

      // check for its size, and bail out if tolerance is achieved
      double sum = 0.0;
      for ( int l=0; l<NumUnknowns; l++ )
          sum += abs(update[l])*abs(update[l]);
      if ( sqrt(sum)<=tol ) {
          converged = true;
          break;
      }

      if (verbose)
          std::cout << "norm_2(psi) = " << sqrt(sum) << std::endl;

      // update
      for ( int l=0; l<NumUnknowns; l++ )
          psi[l] += delta*update[l];

//       if (!k%100) { // add a random component
//           for ( int l=0; l<NumUnknowns; l++ ) {
//               double_complex z = complex<double>( (rand()%100)/1000., (rand()%100)/1000. );
//               psi[l] += z;
//           }
//       }

  }
  // ---------------------------------------------------------------------------

  if (converged)
      std::cout << "Converged!" << std::endl;
  else
      std::cout << "Not converged." << std::endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
