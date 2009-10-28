#include <NOX.H>

#include <Teuchos_DefaultComm.hpp>

#include <EpetraExt_Utils.h> // for toString

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

// User's application specific files
#include "glSystem.h"
#include "ioFactory.h"

// boundary conditions
#include "glBoundaryConditionsInner.h"
#include "glBoundaryConditionsOuter.h"
#include "glBoundaryConditionsCentral.h"

#include <string>

int main(int argc, char *argv[])
{

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Tpetra objects
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm
                                 = Teuchos::DefaultComm<int>::getComm();

  // Get the process ID and the total number of processors
  int MyPID = Comm->getRank();

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

  // ---------------------------------------------------------------------------
  Teuchos::ParameterList problemParameters;
  // define a new dummy psiLexicographic vector, to be adapted instantly
  Teuchos::RCP<Tpetra::Map<int> > dummyMap =
                           Teuchos::rcp ( new Tpetra::Map<int> ( 1, 0, Comm ) );
  Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psiLexicographic =
                                                               Teuchos::ENull();
  if ( withInitialGuess )
    {
      Teuchos::RCP<IoVirtual> fileIo =
               Teuchos::RCP<IoVirtual> ( IoFactory::createFileIo ( filename ) );
      try
        {
          fileIo->read ( psiLexicographic,
			 Comm,
                         problemParameters );
        }
      catch ( const std::exception &e )
        {
          std::cout << e.what() << std::endl;
          return 1;
        }

      if ( psiLexicographic.is_null() ) {
	std::cout << "Input guess empty. Abort." << endl;
	return 1;
      }
    }
  else
    {
      // set the default value
      int    Nx         = 50;
      double edgelength = 10.0;
      double H0         = 0.3;
      std::cout << "Using the standard parameters \n"
                << "    Nx         = " << Nx << ",\n"
                << "    edgelength = " << edgelength << ",\n"
                << "    H0         = " << H0 << "." << std::endl;
      problemParameters.set ( "Nx"        , Nx );
      problemParameters.set ( "edgelength", edgelength );
      problemParameters.set ( "H0"        , H0 );

      int NumGlobalUnknowns = ( Nx+1 ) * ( Nx+1 );
      Teuchos::RCP<Tpetra::Map<int> > standardMap =
           Teuchos::rcp ( new Tpetra::Map<int> ( NumGlobalUnknowns, 0, Comm ) );
      psiLexicographic =
            Teuchos::rcp ( new Tpetra::MultiVector<double_complex,int> ( standardMap,1 ) );
    }
  // ---------------------------------------------------------------------------

  // create the gl problem
  Teuchos::RCP<GlBoundaryConditionsVirtual> boundaryConditions =
		               Teuchos::rcp(new GlBoundaryConditionsCentral() );

  Teuchos::RCP<StaggeredGrid> sGrid =
              Teuchos::rcp( new StaggeredGrid( problemParameters.get<int>("Nx"),
                                               problemParameters.get<double>("edgelength"),
                                               problemParameters.get<double>("H0")
                                             )
                          );

  GinzburgLandau glProblem = GinzburgLandau( sGrid,
                                             boundaryConditions
                                           );

  // ---------------------------------------------------------------------------
  // get proper initial guess
  // ---------------------------------------------------------------------------
  // initialize psi with the same map as psiLexicographic
  Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi =
      Teuchos::rcp ( new Tpetra::MultiVector<double_complex,int> ( psiLexicographic->getMap(),1 ) );
  if ( withInitialGuess )
    {
      // If there was is an initial guess, make sure to get the ordering correct.
      // TODO:
      // Look into having this done by Trilinos. If executed on a multiproc
      // environment, we don't want p to be fully present on all processors.
      int NumComplexUnknowns = glProblem.getStaggeredGrid()->getNumComplexUnknowns();
      std::vector<int> p ( NumComplexUnknowns );
      // fill p:
      glProblem.getStaggeredGrid()->lexicographic2grid ( &p );
      // TODO:
      // The following is certainly not multiproc.
      Teuchos::ArrayRCP<const double_complex> psiView = psiLexicographic->getVector ( 0 )->get1dView();
      for ( int k=0; k<NumComplexUnknowns; k++ )
        {
          psi->replaceGlobalValue ( p[k],
                                    0,
                                    psiView[k]
                                  );
        }
    }
  else
    {
      // create other initial guess for the iteration
      psi->randomize();
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
  double tol       = 1.0e-10;
  bool   converged = false;

  // initialize update value to 0
  bool zeroOut = true;
  Tpetra::MultiVector<double_complex,int> update =
           Tpetra::MultiVector<double_complex,int> ( psi->getMap(),1, zeroOut );
  for( int k=0; k<maxSteps; k++ ) {

      if (verbose) {
          // print the solution to a file
          std::string fileName = outputdir+"/"+"pseudotimestep-"+EpetraExt::toString(k)+".vtk";
          Teuchos::RCP<IoVirtual> fileIo =
               Teuchos::RCP<IoVirtual> ( IoFactory::createFileIo ( fileName ) );
          fileIo->write( *psi,
                          problemParameters,
                         *(glProblem.getStaggeredGrid()) );
      }

      // get the update
      update = glProblem.computeGlVector ( *psi );

      Teuchos::Array<double> norm2array(1);
      update.norm2 ( norm2array() );
      double norm2 = norm2array[0]; 

      // check for its size, and bail out if tolerance is achieved
      if ( norm2<=tol ) {
          converged = true;
          break;
      }

      if (verbose)
          std::cout << "norm_2(psi) = " << norm2 << std::endl;

      // update
      psi->update( 1.0, update, 1.0 );

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
