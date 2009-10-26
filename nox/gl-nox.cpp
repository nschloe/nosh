#include <NOX.H>
#include <NOX_Epetra.H>
// #include <NOX_Abstract_PrePostOperator.H>

#include <Teuchos_DefaultComm.hpp>

// #include <Teuchos_ParameterList.hpp>
// #include <Teuchos_CommandLineProcessor.hpp>

// User's application specific files
#include "ioFactory.h"
#include "ioVtk.h"

#include "EpetraExt_RowMatrixOut.h"

//#include "glPrePostOperator.h"


#include <string>

// for the eigenvalue computation:
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBasicOutputManager.hpp>
#include <AnasaziBlockDavidsonSolMgr.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziLOBPCGSolMgr.hpp>
#include <AnasaziRTRSolMgr.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziConfigDefs.hpp>
#include <AnasaziEpetraAdapter.hpp>
#include <AnasaziBasicSort.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

#include <complex>
typedef std::complex<double> double_complex;

#include "ginzburgLandau.h"
#include "glPrePostOperator.h"
#include "glBoundaryConditionsInner.h"
#include "glBoundaryConditionsOuter.h"
#include "glBoundaryConditionsCentral.h"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

int main ( int argc, char *argv[] )
{

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init ( &argc,&argv );
#endif

  // Create a communicator for Tpetra objects
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm
                                         = Teuchos::DefaultComm<int>::getComm();

  // Get the process ID and the total number of processors
  int MyPID = Comm->getRank();

  // ===========================================================================
  // handle command line arguments
  Teuchos::CommandLineProcessor My_CLP;

  My_CLP.setDocString (
    "This program solves the Ginzburg--Landau problem with a NOX interace.\n"
  );

  bool verbose=false;
  My_CLP.setOption ( "verbose", "silent", &verbose, "Verbostity flag" );

  bool matlabMatrix=false;
  My_CLP.setOption ( "jacobian-file", "no-jacobian-file",&matlabMatrix, "Save the jacobian in a text file" );

  bool computeEigenvalues=false;
  My_CLP.setOption ( "eigenvalues", "no-eigenvalues", &computeEigenvalues, "Compute eigenvalue approximations in the solution" );
  
  bool computeConditionNumber=false;
  My_CLP.setOption ( "condest", "no-condest", &computeConditionNumber, "Compute condition number approximations in the solution" );

  std::string filename = "";
  My_CLP.setOption ( "input-guess", &filename, "File name with initial guess" );

  std::string outputdir = "data";
  My_CLP.setOption ( "output-dir", &outputdir, "Directory to which the solution files are written" );

  // print warning for unrecognized arguments
  My_CLP.recogniseAllOptions ( true );

  // don't throw exceptions
  My_CLP.throwExceptions ( false );

  // finally, parse the stuff!
  Teuchos::CommandLineProcessor::EParseCommandLineReturn
  parseReturn = My_CLP.parse ( argc, argv );
  if ( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED )
    return 0;

  if ( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL )
    return 1; // Error!

  bool withInitialGuess = filename.length() >0;
  // ===========================================================================


  // ---------------------------------------------------------------------------
  Teuchos::ParameterList problemParameters = Teuchos::ParameterList();
  // define a new dummy psiLexicographic vector, to be adapted instantly
  Teuchos::RCP<Tpetra::Map<int> > dummyMap =
    Teuchos::rcp ( new Tpetra::Map<int> ( 1, 0, Comm ) );
  Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psiLexicographic =
                                                               Teuchos::ENull();

//       = Teuchos::rcp( new Tpetra::MultiVector<double_complex,int>(dummyMap,1) );
//       = Teuchos::RCP::;

  if ( withInitialGuess )
    {
      Teuchos::RCP<IoVirtual> fileIo =
               Teuchos::RCP<IoVirtual> ( IoFactory::createFileIo ( filename ) );
      try
        {
          Teuchos::RCP<Teuchos::ParameterList> problemParametersPtr =
                                             Teuchos::rcp( &problemParameters );
          fileIo->read ( psiLexicographic,
			 Comm,
                         problemParametersPtr );
        }
      catch ( const std::exception &e )
        {
          std::cout << e.what() << std::endl;
          return 1;
        }
    }
  else
    {
      // set the default value
      int    Nx         = 50;
      double edgelength = 10.0;
      double H0         = 0.4;
      std::cout << "Using the standard parameters \n"
                << "    Nx         = " << Nx << ",\n"
                << "    edgelength = " << edgelength << ",\n"
                << "    H0         = " << H0 << "." << std::endl;
      problemParameters.set ( "Nx"        , Nx );
      problemParameters.set ( "edgelength", edgelength );
      problemParameters.set ( "H0"        , H0 );

      int NumGlobalUnknowns = ( Nx+1 ) * ( Nx+1 );
      Teuchos::RCP<Tpetra::Map<int> > standardMap
         = Teuchos::rcp ( new Tpetra::Map<int> ( NumGlobalUnknowns, 0, Comm ) );
      psiLexicographic = Teuchos::rcp ( new Tpetra::MultiVector<double_complex,int> ( standardMap,1 ) );
//       psiLexicographic->replaceMap( standardMap );
    }
  // ---------------------------------------------------------------------------

  if ( psiLexicographic.is_null() ) {
    std::cout << "Input guess empty. Abort." << endl;
    return 1; 
  }

  // create the gl problem
  Teuchos::RCP<GlBoundaryConditionsVirtual> boundaryConditions =
                             Teuchos::rcp ( new GlBoundaryConditionsCentral() );

  Teuchos::RCP<StaggeredGrid> sGrid =
                           Teuchos::rcp( new StaggeredGrid( problemParameters.get<int>("Nx"),
                                                            problemParameters.get<double>("edgelength"),
                                                            problemParameters.get<double>("H0")
                                                          )
                                       );

  GinzburgLandau glProblem = GinzburgLandau( sGrid,
                                             boundaryConditions
                                           );


  // create Epetra communicator
#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_MpiComm> eComm
  = Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
  Teuchos::RCP<Epetra_SerialComm>  eComm
  = Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

  // ---------------------------------------------------------------------------
  Teuchos::RCP<GlSystem> glsystem;
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
      Teuchos::RCP<Tpetra::MultiVector<double_complex,int> >  psi
      = Teuchos::rcp ( new Tpetra::MultiVector<double_complex,int> ( psiLexicographic->getMap(),1 ) );
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
      // Create the interface between NOX and the application
      // This object is derived from NOX::Epetra::Interface
      glsystem = Teuchos::rcp ( new GlSystem ( glProblem, eComm, false, psi ) );
    }
  else
    glsystem = Teuchos::rcp ( new GlSystem ( glProblem, eComm, false ) );
  // ---------------------------------------------------------------------------

  // Get initial solution
  Teuchos::RCP<Epetra_Vector> soln = glsystem->getSolution();
  Teuchos::RCP<NOX::Epetra::Vector> noxSoln =
    Teuchos::rcp ( new NOX::Epetra::Vector ( soln,
                   NOX::Epetra::Vector::CreateView ) );


  // ===========================================================================
  // Begin Nonlinear Solver
  // ===========================================================================
  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp ( new Teuchos::ParameterList );
  Teuchos::ParameterList& nlParams = * ( nlParamsPtr.get() );

  // Set the nonlinear solver method
  nlParams.set ( "Nonlinear Solver", "Line Search Based" );

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist ( "Printing" );
  printParams.set ( "MyPID", MyPID );
  printParams.set ( "Output Precision", 16 );
  printParams.set ( "Output Processor", 0 );
  if ( verbose )
    printParams.set ( "Output Information",
                      NOX::Utils::OuterIteration +
                      NOX::Utils::OuterIterationStatusTest +
                      NOX::Utils::InnerIteration +
                      NOX::Utils::LinearSolverDetails +
                      NOX::Utils::Parameters +
                      NOX::Utils::Details +
                      NOX::Utils::Warning +
                      NOX::Utils::Debug +
                      NOX::Utils::TestDetails +
                      NOX::Utils::Error );
  else
    printParams.set ( "Output Information", NOX::Utils::Error +
                      NOX::Utils::TestDetails );

  // Create a print class for controlling output below
  NOX::Utils printing ( printParams );

  // Sublist for line search

  Teuchos::ParameterList& searchParams = nlParams.sublist ( "Line Search" );
  searchParams.set ( "Method", "Full Step" );

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist ( "Direction" );
  dirParams.set ( "Method", "Newton" );
  Teuchos::ParameterList& newtonParams = dirParams.sublist ( "Newton" );
  newtonParams.set ( "Forcing Term Method", "Constant" );

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist ( "Linear Solver" );
//  lsParams.set("Amesos Solver", "Amesos_Superlu");

  // lsParams.set("Aztec Solver", "BiCGStab");
  lsParams.set ( "Aztec Solver", "GMRES" );
  lsParams.set ( "Output Frequency", 32 );
  lsParams.set ( "Output Solver Details", true );
  lsParams.set ( "Max Iterations", 2000 );
  lsParams.set ( "Tolerance", 1e-4 );


  // ====================== BEGIN IFPACK PRECONDITIONER ========================
  lsParams.set ( "Preconditioner","New Ifpack" );
  lsParams.set ( "Ifpack Preconditioner","Amesos" );

  lsParams.set ( "Use Preconditioner as Solver",true );
  lsParams.set ( "Max Age Of Prec", 1 );

  Teuchos::ParameterList& IFPACKparams = lsParams.sublist ( "Ifpack" );
  IFPACKparams.set ( "amesos: solver type","Amesos_Superlu" );
  IFPACKparams.set ( "fact: ilut level-of-fill",1.0 );
  IFPACKparams.set ( "fact: drop tolerance", 1e-1 );
  IFPACKparams.set ( "schwarz: combine mode", "Zero" );
  IFPACKparams.set ( "schwarz: compute condest", true );
  // ====================== END IFPACK PRECONDITIONER ==========================

  // Let's force all status tests to do a full check
  nlParams.sublist ( "Solver Options" ).set ( "Status Test Check Type", "Complete" );

  if ( verbose ) // get custom pre/post actions
    {
      Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo =
        Teuchos::rcp ( new GlPrePostOperator ( glsystem,
                                               problemParameters ) );
      nlParams.sublist ( "Solver Options" )
      .set ( "User Defined Pre/Post Operator", ppo );
    }

  // Create all possible Epetra_Operators.
  Teuchos::RCP<Epetra_RowMatrix> Analytic = glsystem->getJacobian();

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = glsystem;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glsystem;

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
    Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( printParams,
                   lsParams,
                   iReq,
                   iJac,
                   Analytic,
                   *soln ) );

  // Create the Group
  NOX::Epetra::Vector initialGuess ( soln, NOX::Epetra::Vector::CreateView );
  Teuchos::RCP<NOX::Epetra::Group> grpPtr =
    Teuchos::rcp ( new NOX::Epetra::Group ( printParams,
                                            iReq,
                                            initialGuess,
                                            linSys ) );

  NOX::Epetra::Group& grp = *grpPtr;
  // ------------------------------------------------------------------------


  // ------------------------------------------------------------------------
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp ( new NOX::StatusTest::NormF ( 1.0e-10 ) );
  Teuchos::RCP<NOX::StatusTest::NormF> relresid =
    Teuchos::rcp ( new NOX::StatusTest::NormF ( grp, 1.0e-10 ) );
  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp ( new NOX::StatusTest::NormUpdate ( 1.0e-10 ) );
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp ( new NOX::StatusTest::NormWRMS ( 1.0e-2, 1.0e-5 ) );
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp ( new NOX::StatusTest::Combo ( NOX::StatusTest::Combo::AND ) );
  converged->addStatusTest ( absresid );
  converged->addStatusTest ( relresid );
  converged->addStatusTest ( wrms );
  converged->addStatusTest ( update );

  int maxNonlinearIterations = 500;
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp ( new NOX::StatusTest::MaxIters ( maxNonlinearIterations ) );
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp ( new NOX::StatusTest::FiniteValue );
  // this test is useful if we start in a solution
  Teuchos::RCP<NOX::StatusTest::NormF> absresexact =
    Teuchos::rcp ( new NOX::StatusTest::NormF ( 1.0e-13 ) );
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp ( new NOX::StatusTest::Combo ( NOX::StatusTest::Combo::OR ) );

  combo->addStatusTest ( fv );
  combo->addStatusTest ( absresexact );
  combo->addStatusTest ( converged );
  combo->addStatusTest ( maxiters );
  // ------------------------------------------------------------------------

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver ( grpPtr,
      combo,
      nlParamsPtr );
  // solve!
  NOX::StatusTest::StatusType solvStatus = solver->solve();
  // ===========================================================================
  // End Nonlinear Solver
  // ===========================================================================

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&> ( solver->getSolutionGroup() );
  const Epetra_Vector& finalSolution =
    ( dynamic_cast<const NOX::Epetra::Vector&> ( finalGroup.getX() ) ).
    getEpetraVector();

  if ( computeConditionNumber )
    /*        {
      // -----------------------------------------------------------------------
      // compute the condition number
      try
        {
          grpPtr->computeJacobian();
          grpPtr->computeJacobianConditionNumber ( 2000, 1e-2, 30, true );
          double kappa = finalGroup.getJacobianConditionNumber();
          std::cout << "Condition number: kappa = " << kappa << "." << std::endl;
        }
      catch ( std::exception& e )
        {
          std::cerr << e.what() << std::endl;
        }
      // -----------------------------------------------------------------------
      }*/
  if ( matlabMatrix )
    {
      std::string jacFilename = "jacobianMatrix.dat";
      EpetraExt::RowMatrixToMatlabFile(jacFilename.c_str(),*(glsystem->getJacobian()));
    }

  if ( computeEigenvalues )
    {
      // -----------------------------------------------------------------------
      // compute some eigenvalues using anasazi
      bool debug = true; // even more increased verbosity
      std::string which ( "LR" ); // compute the rightmost eigenvalues

      bool boolret;

      typedef double ScalarType;
      typedef Teuchos::ScalarTraits<ScalarType>          SCT;
      typedef SCT::magnitudeType               MagnitudeType;
      typedef Epetra_MultiVector                          MV;
      typedef Epetra_Operator                             OP;
      typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
      typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;

      Epetra_BlockMap Map = finalSolution.Map();

      // shortname for the final Jacobian
      Teuchos::RCP<Epetra_Operator> J = grpPtr->getLinearSystem()->getJacobianOperator();


      //Teuchos::RCP<Epetra_CrsMatrix> Jcrs = grpPtr->getLinearSystem()->getJacobianOperator();


      // - - - - - - - - - - - - - - - - -
      // Start the block Arnoldi iteration
      // - - - - - - - - - - - - - - - - -
      // Variables used for the Block Krylov Schur Method
      int nev = 10;
      int blockSize = 1;
      int numBlocks = 30;
      int maxRestarts = 250;
      double tol = 1e-10;

      // Create a sort manager to pass into the block Krylov-Schur solver manager
      // -->  Make sure the reference-counted pointer is of type Anasazi::SortManager<>
      // -->  The block Krylov-Schur solver manager uses Anasazi::BasicSort<> by default,
      //      so you can also pass in the parameter "Which", instead of a sort manager.
      Teuchos::RCP<Anasazi::SortManager<MagnitudeType> > MySort =
               Teuchos::rcp ( new Anasazi::BasicSort<MagnitudeType> ( which ) );

      // Set verbosity level
      int verbosity = Anasazi::Errors + Anasazi::Warnings;
      verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
      if ( debug )
        {
          verbosity += Anasazi::Debug;
        }

      // Create parameter list to pass into solver manager
      Teuchos::ParameterList MyPL;
      MyPL.set ( "Verbosity", verbosity );
      MyPL.set ( "Sort Manager", MySort );
      //MyPL.set( "Which", which );
      MyPL.set ( "Block Size", blockSize );
      MyPL.set ( "Num Blocks", numBlocks );
      MyPL.set ( "Maximum Restarts", maxRestarts );
      //MyPL.set( "Step Size", stepSize );
      MyPL.set ( "Convergence Tolerance", tol );

      // Create an Epetra_MultiVector for an initial vector to start the solver.
      // Note:  This needs to have the same number of columns as the blocksize.
      Teuchos::RCP<Epetra_MultiVector> ivec =
                     Teuchos::rcp ( new Epetra_MultiVector ( Map, blockSize ) );
      ivec->Random();

      // Create the eigenproblem.
      Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
        Teuchos::rcp ( new Anasazi::BasicEigenproblem<double, MV, OP> ( J, ivec ) );

      // Set the number of eigenvalues requested
      MyProblem->setNEV ( nev );

      // Inform the eigenproblem that you are finishing passing it information
      boolret = MyProblem->setProblem();
      if ( boolret != true )
        {
          if ( MyPID == 0 )
            {
              cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
            }
#ifdef HAVE_MPI
          MPI_Finalize() ;
#endif
          return -1;
        }

      // Initialize the Block Arnoldi solver
//       Anasazi::BlockDavidsonSolMgr<double, MV, OP> 
//       Anasazi::LOBPCGSolMgr<double, MV, OP> 
//       Anasazi::RTRSolMgr<double, MV, OP> 
      Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> 
                                                MySolverMgr ( MyProblem, MyPL );

      // Solve the problem to the specified tolerances or length
      Anasazi::ReturnType returnCode = MySolverMgr.solve();
      if ( returnCode != Anasazi::Converged && MyPID==0 && verbose )
        cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;

      // Get the Ritz values from the eigensolver
//       std::vector<Anasazi::Value<double> > ritzValues = MySolverMgr.getRitzValues();

//       // Output computed eigenvalues and their direct residuals
//       if (MyPID==0) {
//         int numritz = (int)ritzValues.size();
//         cout.setf(std::ios_base::right, std::ios_base::adjustfield);
//         cout<<endl<< "Computed Ritz Values"<< endl;
//         if (MyProblem->isHermitian()) {
//           cout<< std::setw(16) << "Real Part"
//             << endl;
//           cout<<"-----------------------------------------------------------"<<endl;
//           for (int i=0; i<numritz; i++) {
//             cout<< std::setw(16) << ritzValues[i].realpart
//               << endl;
//           }
//           cout<<"-----------------------------------------------------------"<<endl;
//         }
//         else {
//           cout<< std::setw(16) << "Real Part"
//             << std::setw(16) << "Imag Part"
//             << endl;
//           cout<<"-----------------------------------------------------------"<<endl;
//           for (int i=0; i<numritz; i++) {
//             cout<< std::setw(16) << ritzValues[i].realpart
//               << std::setw(16) << ritzValues[i].imagpart
//               << endl;
//           }
//           cout<<"-----------------------------------------------------------"<<endl;
//           }
//         }

      // Get the eigenvalues and eigenvectors from the eigenproblem
      Anasazi::Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
      std::vector<Anasazi::Value<ScalarType> > evals = sol.Evals;
      Teuchos::RCP<MV> evecs = sol.Evecs;
      std::vector<int> index = sol.index;
      int numev = sol.numVecs;

      Teuchos::RCP<Epetra_Vector> ev1 ( ( *evecs ) ( 0 ) );

//   glsystem->solutionToFile( *ev1,
//                             problemParameters,
//                             "data/ev1.vtk" );
//
//       cout<< std::setw(16) << "Real Part"
//           << std::setw(16) << "Imag Part" << endl;
//       cout<<"-----------------------------------------------------------"<<endl;
//       for (int i=0; i<numev; i++) {
//           cout<< std::setw(16) << evals[i].realpart
//             << std::setw(16) << evals[i].imagpart << endl;
//       }
      // -----------------------------------------------------------------------
    }

  // ---------------------------------------------------------------------------
  // print the solution to a file
  //problemParameters.set( "FE", glsystem.freeEnergy
  std::string solutionFile = "data/solution.vtk";
  glsystem->solutionToFile ( finalSolution,
                             problemParameters,
                             solutionFile );
  // ---------------------------------------------------------------------------


  // Tests
  int status = 0; // Converged

  // 1. Convergence
  if ( solvStatus != NOX::StatusTest::Converged )
    {
      status = 1;
      if ( printing.isPrintType ( NOX::Utils::Error ) )
        printing.out() << "Nonlinear solver failed to converge!" << endl;
    }
#ifndef HAVE_MPI
  // 2. Linear solve iterations (53) - SERIAL TEST ONLY!
  //    The number of linear iterations changes with # of procs.
  if ( const_cast<Teuchos::ParameterList&> ( solver->getList() ).sublist ( "Direction" ).sublist ( "Newton" ).sublist ( "Linear Solver" ).sublist ( "Output" ).get ( "Total Number of Linear Iterations",0 ) != 53 )
    {
      status = 2;
    }
#endif
  // 3. Nonlinear solve iterations (10)
  if ( const_cast<Teuchos::ParameterList&> ( solver->getList() ).sublist ( "Output" ).get ( "Nonlinear Iterations", 0 ) == maxNonlinearIterations )
    status = 3;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // Final return value (0 = successful, non-zero = failure)
  return status;
}
