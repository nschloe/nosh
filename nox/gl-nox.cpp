#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Abstract_PrePostOperator.H>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

// User's application specific files
#include "glSystem.h"
#include "ioFactory.h"

#include "glPrePostOperator.h"

#include <string>

// for the eigenvalue computation:
// #include <AnasaziBasicEigenproblem.hpp>
// #include <AnasaziBasicOutputManager.hpp>
// #include <AnasaziBlockDavidsonSolMgr.hpp>

#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicSort.hpp"

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

  if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL )
    return 1; // Error!

  bool withInitialGuess = filename.length()>0;
  // ===========================================================================


  // ---------------------------------------------------------------------------
  std::vector<double_complex> psiLexicographic;
  Teuchos::ParameterList      problemParameters;
  if (withInitialGuess) {
      Teuchos::RCP<IoVirtual> fileIo = Teuchos::RCP<IoVirtual>( IoFactory::createFileIo( filename ) );
      try {
          fileIo->read( &psiLexicographic,
                        &problemParameters );
      }
      catch ( const std::exception &e ) {
          std::cout << e.what() << std::endl;
          return 1;
      }
  } else {
      // set the default value
      int Nx = 10;
      double edgelength = 10.0;
      double H0 = 0.4;
      std::cout << "Using the standard parameters \n"
                << "    Nx         = " << Nx << ",\n"
                << "    edgelength = " << edgelength << ",\n"
                << "    H0         = " << H0 << "." << std::endl;
      problemParameters.set( "Nx"        , Nx   );
      problemParameters.set( "edgelength", edgelength );
      problemParameters.set( "H0"        , H0  );
  }
  // ---------------------------------------------------------------------------

  // create the gl problem
  GinzburgLandau glProblem = GinzburgLandau( problemParameters.get<int>("Nx"),
                                             problemParameters.get<double>("edgelength"),
                                             problemParameters.get<double>("H0")
                                           );

  // ---------------------------------------------------------------------------
  Teuchos::RCP<GlSystem> glsystem;
  if (withInitialGuess) {
      // If there was is an initial guess, make sure to get the ordering correct.
      int              NumUnknowns = glProblem. getStaggeredGrid()
                                              ->getNumComplexUnknowns();
      std::vector<int> p(NumUnknowns);
      // fill p:
      glProblem.getStaggeredGrid()->lexicographic2grid( &p );
      std::vector<double_complex> psi(NumUnknowns);
      for (int k=0; k<NumUnknowns; k++)
          psi[p[k]] = psiLexicographic[k];

      // Create the interface between NOX and the application
      // This object is derived from NOX::Epetra::Interface
      glsystem = Teuchos::rcp(new GlSystem( glProblem, Comm, &psi ) );
  } else
      glsystem = Teuchos::rcp(new GlSystem( glProblem,Comm ) );
  // ---------------------------------------------------------------------------


  // Get initial solution
  Teuchos::RCP<Epetra_Vector> soln = glsystem->getSolution();
  Teuchos::RCP<NOX::Epetra::Vector> noxSoln =
         Teuchos::rcp(new NOX::Epetra::Vector(soln,
                                              NOX::Epetra::Vector::CreateView));


  // ===========================================================================
  // Begin Nonlinear Solver
  // ===========================================================================
  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID);
  printParams.set("Output Precision", 16);
  printParams.set("Output Processor", 0);
  if (verbose)
    printParams.set("Output Information",
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
    printParams.set("Output Information", NOX::Utils::Error +
                                          NOX::Utils::TestDetails );

  // Create a print class for controlling output below
  NOX::Utils printing(printParams);

  // Sublist for line search

  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
//  lsParams.set("Amesos Solver", "Amesos_Superlu");

  // lsParams.set("Aztec Solver", "BiCGStab");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Output Frequency", 32);
  lsParams.set("Output Solver Details", true);
  lsParams.set("Max Iterations", 2000);
  lsParams.set("Tolerance", 1e-4);


  // ====================== BEGIN IFPACK PRECONDITIONER ========================
  lsParams.set("Preconditioner","New Ifpack");
  lsParams.set("Ifpack Preconditioner","Amesos");

  lsParams.set("Use Preconditioner as Solver",true);
  lsParams.set("Max Age Of Prec", 1);

  Teuchos::ParameterList& IFPACKparams = lsParams.sublist("Ifpack");
      IFPACKparams.set("amesos: solver type","Amesos_Superlu");
      IFPACKparams.set("fact: ilut level-of-fill",1.0);
      IFPACKparams.set("fact: drop tolerance", 1e-1);
      IFPACKparams.set("schwarz: combine mode", "Zero");
      IFPACKparams.set("schwarz: compute condest", true);
  // ====================== END IFPACK PRECONDITIONER ==========================

  // Let's force all status tests to do a full check
  nlParams.sublist("Solver Options").set("Status Test Check Type", "Complete");

  if (verbose) { // get custom pre/post actions
      Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo =
                                      Teuchos::rcp(new GlPrePostOperator(glsystem,
                                                                         problemParameters));
      nlParams.sublist("Solver Options")
              .set("User Defined Pre/Post Operator", ppo);
  }

  // Create all possible Epetra_Operators.
  Teuchos::RCP<Epetra_RowMatrix> Analytic = glsystem->getJacobian();

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = glsystem;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glsystem;

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
           Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO( printParams,
                                                              lsParams,
                                                              iReq,
                                                              iJac,
                                                              Analytic,
                                                              *soln        ) );

  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::Epetra::Vector::CreateView);
  Teuchos::RCP<NOX::Epetra::Group> grpPtr =
    Teuchos::rcp(new NOX::Epetra::Group(printParams,
                                        iReq,
                                        initialGuess,
                                        linSys));

  NOX::Epetra::Group& grp = *grpPtr;
  // ------------------------------------------------------------------------


  // ------------------------------------------------------------------------
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-10));
  Teuchos::RCP<NOX::StatusTest::NormF> relresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-10));
  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-10));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-5));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(relresid);
  converged->addStatusTest(wrms);
  converged->addStatusTest(update);

  int maxNonlinearIterations = 500;
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNonlinearIterations));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  // this test is useful if we start in a solution
  Teuchos::RCP<NOX::StatusTest::NormF> absresexact =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-14));
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  combo->addStatusTest(fv);
  combo->addStatusTest(absresexact);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);
  // ------------------------------------------------------------------------

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver( grpPtr,
                                                                        combo,
                                                                        nlParamsPtr );
  // solve!
  NOX::StatusTest::StatusType solvStatus = solver->solve();
  // ===========================================================================
  // End Nonlinear Solver
  // ===========================================================================

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& finalSolution =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).
    getEpetraVector();

  if (verbose) {
      try {
          // get condition number
          grpPtr->computeJacobian();
          grpPtr->computeJacobianConditionNumber( 2000, 1e-2, 30, true );
          double kappa = finalGroup.getJacobianConditionNumber();
          std::cout << "Condition number: kappa = " << kappa << "." << std::endl;
      }
      catch ( std::exception& e ) {
          std::cerr << e.what() << std::endl;
      }

//       // check back with
//       NOX::Epetra::Vector nullVec(finalSolution);
//       double tmp;
//       // Construct the (complex) vector i*finalSolution.
//       int k=0;
//       Epetra_Vector & nV = nullVec.getEpetraVector();
//       while ( k<nullVec.length() ) {
//           tmp     =  nV[k];
//           nV[k]   = -nV[k+1];
//           nV[k+1] = tmp;
//           k += 2;
//       }
//       NOX::Epetra::Vector resVec(finalSolution);
//       grpPtr->applyJacobian( nullVec, resVec );
//       double norm;
//       resVec.getEpetraVector().Norm2( &norm );
//       std::cout << "Norm: ||v|| = " << norm << std::endl;

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  bool verbose = true;
  bool debug = false;
  std::string which("SM");

  bool boolret;

  typedef double ScalarType;
  typedef Teuchos::ScalarTraits<ScalarType>          SCT;
  typedef SCT::magnitudeType               MagnitudeType;
  typedef Epetra_MultiVector                          MV;
  typedef Epetra_Operator                             OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;

  Epetra_BlockMap Map = finalSolution.Map();

  Teuchos::RCP<Epetra_Operator> J = grpPtr->getLinearSystem()->getJacobianOperator();

  //************************************
  // Start the block Arnoldi iteration
  //***********************************
  //
  //  Variables used for the Block Krylov Schur Method
  //    
  int nev = 4;
  int blockSize = 1;
  int numBlocks = 20;
  int maxRestarts = 100;
  //int stepSize = 5;
  double tol = 1e-8;

  // Create a sort manager to pass into the block Krylov-Schur solver manager
  // -->  Make sure the reference-counted pointer is of type Anasazi::SortManager<>
  // -->  The block Krylov-Schur solver manager uses Anasazi::BasicSort<> by default,
  //      so you can also pass in the parameter "Which", instead of a sort manager.
  Teuchos::RCP<Anasazi::SortManager<MagnitudeType> > MySort =     
    Teuchos::rcp( new Anasazi::BasicSort<MagnitudeType>( which ) );

  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }
  //
  // Create parameter list to pass into solver manager
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Sort Manager", MySort );
  //MyPL.set( "Which", which );  
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  //MyPL.set( "Step Size", stepSize );
  MyPL.set( "Convergence Tolerance", tol );

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();

  // Create the eigenproblem.
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(J, ivec) );
  
//   // Inform the eigenproblem that the operator A is symmetric
//   MyProblem->setHermitian(rho==0.0); 
  
  // Set the number of eigenvalues requested
  MyProblem->setNEV( nev );
  
  // Inform the eigenproblem that you are finishing passing it information
  boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }
  
  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);
  
  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && MyPID==0 && verbose) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
  }

//   // Get the Ritz values from the eigensolver
//   std::vector<Anasazi::Value<double> > ritzValues = MySolverMgr.getRitzValues();
//   
//   // Output computed eigenvalues and their direct residuals
//   if (verbose && MyPID==0) {
//     int numritz = (int)ritzValues.size();
//     cout.setf(std::ios_base::right, std::ios_base::adjustfield);
//     cout<<endl<< "Computed Ritz Values"<< endl;
//     if (MyProblem->isHermitian()) {
//       cout<< std::setw(16) << "Real Part"
//         << endl;
//       cout<<"-----------------------------------------------------------"<<endl;
//       for (int i=0; i<numritz; i++) {
//         cout<< std::setw(16) << ritzValues[i].realpart 
//           << endl;
//       }  
//       cout<<"-----------------------------------------------------------"<<endl;
//     } 
//     else {
//       cout<< std::setw(16) << "Real Part"
//         << std::setw(16) << "Imag Part"
//         << endl;
//       cout<<"-----------------------------------------------------------"<<endl;
//       for (int i=0; i<numritz; i++) {
//         cout<< std::setw(16) << ritzValues[i].realpart 
//           << std::setw(16) << ritzValues[i].imagpart 
//           << endl;
//       }  
//       cout<<"-----------------------------------------------------------"<<endl;
//       }  
//     }

//   // Get the eigenvalues and eigenvectors from the eigenproblem
//   Anasazi::Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
//   std::vector<Anasazi::Value<ScalarType> > evals = sol.Evals;
//   Teuchos::RCP<MV> evecs = sol.Evecs;
//   std::vector<int> index = sol.index;
//   int numev = sol.numVecs;
//   
//   if (numev > 0) {
//     // Compute residuals.
//     Teuchos::LAPACK<int,double> lapack;
//     std::vector<double> normA(numev);
//     
//     if (MyProblem->isHermitian()) {
//       // Get storage
//       Epetra_MultiVector Aevecs(Map,numev);
//       Teuchos::SerialDenseMatrix<int,double> B(numev,numev);
//       B.putScalar(0.0); 
//       for (int i=0; i<numev; i++) {B(i,i) = evals[i].realpart;}
//       
//       // Compute A*evecs
//       OPT::Apply( *J, *evecs, Aevecs );
//       
//       // Compute A*evecs - lambda*evecs and its norm
//       MVT::MvTimesMatAddMv( -1.0, *evecs, B, 1.0, Aevecs );
//       MVT::MvNorm( Aevecs, normA );
//       
//       // Scale the norms by the eigenvalue
//       for (int i=0; i<numev; i++) {
//         normA[i] /= Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
//       }
//     } else {
//       // The problem is non-Hermitian.
//       int i=0;
//       std::vector<int> curind(1);
//       std::vector<double> resnorm(1), tempnrm(1);
//       Teuchos::RCP<MV> evecr, eveci, tempAevec;
//       Epetra_MultiVector Aevec(Map,numev);
//       
//       // Compute A*evecs
//       OPT::Apply( *J, *evecs, Aevec );
//       
//       Teuchos::SerialDenseMatrix<int,double> Breal(1,1), Bimag(1,1);
//       while (i<numev) {
//         if (index[i]==0) {
//           // Get a view of the current eigenvector (evecr)
//           curind[0] = i;
//           evecr = MVT::CloneView( *evecs, curind );
// 
//           // Get a copy of A*evecr
//           tempAevec = MVT::CloneCopy( Aevec, curind );
// 
//           // Compute A*evecr - lambda*evecr
//           Breal(0,0) = evals[i].realpart;
//           MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
// 
//           // Compute the norm of the residual and increment counter
//           MVT::MvNorm( *tempAevec, resnorm );
//           normA[i] = resnorm[0]/Teuchos::ScalarTraits<MagnitudeType>::magnitude( evals[i].realpart );
//           i++;
//         } else {
//           // Get a view of the real part of the eigenvector (evecr)
//           curind[0] = i;
//           evecr = MVT::CloneView( *evecs, curind );
// 
//           // Get a copy of A*evecr
//           tempAevec = MVT::CloneCopy( Aevec, curind );
// 
//           // Get a view of the imaginary part of the eigenvector (eveci)
//           curind[0] = i+1;
//           eveci = MVT::CloneView( *evecs, curind );
// 
//           // Set the eigenvalue into Breal and Bimag
//           Breal(0,0) = evals[i].realpart;
//           Bimag(0,0) = evals[i].imagpart;
// 
//           // Compute A*evecr - evecr*lambdar + eveci*lambdai
//           MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
//           MVT::MvTimesMatAddMv( 1.0, *eveci, Bimag, 1.0, *tempAevec );
//           MVT::MvNorm( *tempAevec, tempnrm );
// 
//           // Get a copy of A*eveci
//           tempAevec = MVT::CloneCopy( Aevec, curind );
// 
//           // Compute A*eveci - eveci*lambdar - evecr*lambdai
//           MVT::MvTimesMatAddMv( -1.0, *evecr, Bimag, 1.0, *tempAevec );
//           MVT::MvTimesMatAddMv( -1.0, *eveci, Breal, 1.0, *tempAevec );
//           MVT::MvNorm( *tempAevec, resnorm );
// 
//           // Compute the norms and scale by magnitude of eigenvalue
//           normA[i] = lapack.LAPY2( tempnrm[i], resnorm[i] ) /
//             lapack.LAPY2( evals[i].realpart, evals[i].imagpart );
//           normA[i+1] = normA[i];
// 
//           i=i+2;
//         }
//       }
//     }
// 
//     // Output computed eigenvalues and their direct residuals
//     if (verbose && MyPID==0) {
//       cout.setf(std::ios_base::right, std::ios_base::adjustfield);
//       cout<<endl<< "Actual Residuals"<<endl;
//       if (MyProblem->isHermitian()) {
//         cout<< std::setw(16) << "Real Part"
//           << std::setw(20) << "Direct Residual"<< endl;
//         cout<<"-----------------------------------------------------------"<<endl;
//         for (int i=0; i<numev; i++) {
//           cout<< std::setw(16) << evals[i].realpart 
//             << std::setw(20) << normA[i] << endl;
//         }  
//         cout<<"-----------------------------------------------------------"<<endl;
//       } 
//       else {
//         cout<< std::setw(16) << "Real Part"
//           << std::setw(16) << "Imag Part"
//           << std::setw(20) << "Direct Residual"<< endl;
//         cout<<"-----------------------------------------------------------"<<endl;
//         for (int i=0; i<numev; i++) {
//           cout<< std::setw(16) << evals[i].realpart 
//             << std::setw(16) << evals[i].imagpart 
//             << std::setw(20) << normA[i] << endl;
//         }  
//         cout<<"-----------------------------------------------------------"<<endl;
//       }  
//     }
//   }
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ---------------------------------------------------------------------------
  // get eigenvalue with largest real part with Anasazi

//   // Create an Anasazi output manager
//   Anasazi::BasicOutputManager<double> printer;
//   printer.stream(Anasazi::Errors) << Anasazi::Anasazi_Version() << std::endl << std::endl;
// 
//   int blockSize = finalSolution.GlobalLength();
//   int numBlocks = 1;
// 
//   // Create an Epetra_MultiVector for an initial vector to start the solver.
//   // Note:  This needs to have the same number of columns as the blocksize.
//   Epetra_BlockMap vecMap = finalSolution.Map();
//   Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(vecMap, blockSize) );
//   // initialize as random
//   ivec->Random();
// 
//   // get pointer to Jacobian
//   Teuchos::RCP<Epetra_Operator> J = grpPtr->getLinearSystem()->getJacobianOperator();
// 
//   // Create the eigenproblem.
//   Teuchos::RCP<Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator> > MyProblem =
//     Teuchos::rcp( new Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator>(J,ivec) );
// 
//   // Set the number of eigenvalues requested
//   //
//   int nev = 4;
//   MyProblem->setNEV( nev );
// 
//   // Inform the eigenproblem that you are finishing passing it information
//   bool boolret = MyProblem->setProblem();
//   if (boolret != true) {
//     printer.print(Anasazi::Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
//   }
// 
//   // Create parameter list to pass into the solver manager
//   Teuchos::ParameterList MyPL;
//   std::string which("LR");  // get eigenvalues with the largest real part
//   MyPL.set( "Which", which );
//   MyPL.set( "Block Size", blockSize );
//   MyPL.set( "Num Blocks", numBlocks );
//   int maxRestarts = 100;
//   MyPL.set( "Maximum Restarts", maxRestarts );
//   double tol = 1e-6;
//   MyPL.set( "Convergence Tolerance", tol );
// 
//   // create the solver manager
//   Anasazi::BlockDavidsonSolMgr<double, Epetra_MultiVector, Epetra_Operator> MySolverMan(MyProblem, MyPL);

//   // Solve the problem
//   Anasazi::ReturnType returnCode = MySolverMan.solve();
// 
//   // Get the eigenvalues and eigenvectors from the eigenproblem
//   Anasazi::Eigensolution<double,Epetra_MultiVector> sol = MyProblem->getSolution();
//   std::vector<Anasazi::Value<double> > evals = sol.Evals;
//   Teuchos::RCP<Epetra_MultiVector> evecs = sol.Evecs;
// 
//   // Compute residuals.
//   std::vector<double> normR(sol.numVecs);
//   if (sol.numVecs > 0) {
//     Teuchos::SerialDenseMatrix<int,double> T(sol.numVecs, sol.numVecs);
//     Epetra_MultiVector tempAevec( vecMap, sol.numVecs );
//     T.putScalar(0.0);
//     for (int i=0; i<sol.numVecs; i++) {
//       T(i,i) = evals[i].realpart;
//     }
//     J->Apply( *evecs, tempAevec );
//     MVT::MvTimesMatAddMv( -1.0, *evecs, T, 1.0, tempAevec );
//     MVT::MvNorm( tempAevec, normR );
//   }
// 
//   // Print the results
//   std::ostringstream os;
//   os.setf(std::ios_base::right, std::ios_base::adjustfield);
//   os<<"Solver manager returned "
//     << (returnCode == Anasazi::Converged ? "converged." : "unconverged.")
//     << std::endl;
//   os<<std::endl;
//   os<<"------------------------------------------------------"<<std::endl;
//   os<<std::setw(16)<<"Eigenvalue"
//     <<std::setw(18)<<"Direct Residual"
//     <<std::endl;
//   os<<"------------------------------------------------------"<<std::endl;
//   for (int i=0; i<sol.numVecs; i++) {
//     os<<std::setw(16)<<evals[i].realpart
//       <<std::setw(18)<<normR[i]/evals[i].realpart
//       <<std::endl;
//   }
//   os<<"------------------------------------------------------"<<std::endl;
//   printer.print(Anasazi::Errors,os.str());
//   // ---------------------------------------------------------------------------
  }

  // ---------------------------------------------------------------------------
  // print the solution to a file
  glsystem->solutionToFile( finalSolution,
                            problemParameters,
                            "data/solution.vtk" );
  // ---------------------------------------------------------------------------


  // Tests
  int status = 0; // Converged

  // 1. Convergence
  if (solvStatus != NOX::StatusTest::Converged) {
      status = 1;
      if (printing.isPrintType(NOX::Utils::Error))
        printing.out() << "Nonlinear solver failed to converge!" << endl;
  }
#ifndef HAVE_MPI
  // 2. Linear solve iterations (53) - SERIAL TEST ONLY!
  //    The number of linear iterations changes with # of procs.
  if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").get("Total Number of Linear Iterations",0) != 53) {
    status = 2;
  }
#endif
  // 3. Nonlinear solve iterations (10)
  if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Output").get("Nonlinear Iterations", 0) == maxNonlinearIterations)
    status = 3;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // Final return value (0 = successful, non-zero = failure)
  return status;
}
