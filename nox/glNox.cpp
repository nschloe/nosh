#include "glNox.h"

#include "ioVirtual.h"
#include "ioFactory.h"
#include "glException.h"
#include "glBoundaryConditionsVirtual.h"
// #include "glBoundaryConditionsInner.h"
// #include "glBoundaryConditionsOuter.h"
#include "glBoundaryConditionsCentral.h"
#include "ginzburgLandau.h"
#include "glPrePostOperator.h"

#include <NOX.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_Epetra.H>

// // #include <NOX_Abstract_PrePostOperator.H>
// 
// // #include <Teuchos_ParameterList.hpp>
// 
// // User's application specific files
// #include "ioVtk.h"
// 
// //#include "glPrePostOperator.h"
// 
// #include <Teuchos_RCP.hpp>
// #include <Teuchos_CommandLineProcessor.hpp>
// #include <Teuchos_ParameterList.hpp>
// 
// #include <Tpetra_Map.hpp>
// #include <Tpetra_MultiVector.hpp>
// 
// #include <complex>
// typedef std::complex<double> double_complex;
// 
// #include "glPrePostOperator.h"


// for the eigenvalue computation:
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBasicOutputManager.hpp>
#include <AnasaziBlockDavidsonSolMgr.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziLOBPCGSolMgr.hpp>
#include <AnasaziRTRSolMgr.hpp>
#include <AnasaziConfigDefs.hpp>
#include <AnasaziEpetraAdapter.hpp>
#include <AnasaziBasicSort.hpp>

// =============================================================================
glNox::glNox( const std::string fileName,
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
              const Teuchos::RCP<const Epetra_Comm>         &eComm ) :
  Comm_( comm ),
  eComm_( eComm ),
  problemParameters_(),
  MyPID_( comm->getRank() ),
  glSystem_( Teuchos::ENull() ),
  nlParamsPtr_( Teuchos::rcp ( new Teuchos::ParameterList ) ),
  combo_( Teuchos::ENull() ),
  solver_( Teuchos::ENull() ),
  verbose_( false ),
  maxNonlinearIterations_( 10 )
{
  // instantiate file I/O object
  Teuchos::RCP<IoVirtual> fileIo =
	       Teuchos::RCP<IoVirtual> ( IoFactory::createFileIo ( fileName ) );

  Teuchos::RCP<Tpetra::Vector<double_complex,int> > psi = Teuchos::ENull();

  // read the stuff
  fileIo->read ( psi,
                 Comm_,
                 problemParameters_ );

  if ( psi.is_null() )
    throw glException( "glNox::glNox", "Input guess empty" );

  Teuchos::RCP<GlBoundaryConditionsVirtual> boundaryConditions =
                             Teuchos::rcp ( new GlBoundaryConditionsCentral() );

  int    Nx         = problemParameters_.get<int>   ("Nx");
  double edgeLength = problemParameters_.get<double>("edgelength");
  double H0         = problemParameters_.get<double>("H0");
  Teuchos::RCP<StaggeredGrid> sGrid =
                       Teuchos::rcp ( new StaggeredGrid( Nx, edgeLength, H0 ) );  

  GinzburgLandau glProblem = GinzburgLandau( sGrid,
                                             boundaryConditions
                                           );

  reOrder( *psi, sGrid );

  // Create the interface between NOX and the application
  // This object is derived from NOX::Epetra::Interface
  bool reverse = false;
  glSystem_ = Teuchos::rcp ( new GlSystem ( glProblem, eComm, reverse, psi ) );
}
// =============================================================================
glNox::glNox( const int Nx,
              const double edgeLength,
              const double H0,
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
              const Teuchos::RCP<const Epetra_Comm>         &eComm ) :
  Comm_( comm ),
  eComm_( eComm ),
  problemParameters_(),
  MyPID_( comm->getRank() ),
  glSystem_( Teuchos::ENull() ),
  nlParamsPtr_( Teuchos::rcp ( new Teuchos::ParameterList ) ),
  combo_( Teuchos::ENull() ),
  solver_( Teuchos::ENull() ),
  verbose_( false ),
  maxNonlinearIterations_( 10 )
{
  problemParameters_.set ( "Nx"        , Nx );
  problemParameters_.set ( "edgelength", edgeLength );
  problemParameters_.set ( "H0"        , H0 );

  Teuchos::RCP<GlBoundaryConditionsVirtual> boundaryConditions =
                             Teuchos::rcp ( new GlBoundaryConditionsCentral() );

  Teuchos::RCP<StaggeredGrid> sGrid =
                       Teuchos::rcp ( new StaggeredGrid( Nx, edgeLength, H0 ) );  

  GinzburgLandau glProblem = GinzburgLandau( sGrid,
                                             boundaryConditions
                                           );

  bool reverse = false;
  glSystem_ = Teuchos::rcp ( new GlSystem ( glProblem, eComm, reverse ) );
}
// =============================================================================
void
glNox::solve()
{
  NOX::StatusTest::StatusType solvStatus = solver_->solve();
}
// =============================================================================
// set parameters
void
glNox::setSolverOptions( int maxNonlinearIterations )
{
  maxNonlinearIterations_ = maxNonlinearIterations;

  nlParamsPtr_ =  Teuchos::rcp ( new Teuchos::ParameterList );

  Teuchos::ParameterList& nlParams = * ( nlParamsPtr_.get() );
  setNonlinearSolverParameters( nlParams );
}
// =============================================================================
void
glNox::createSolverGroup()
{
  // Create all possible Epetra_Operators.
  Teuchos::RCP<Epetra_RowMatrix> Analytic = glSystem_->getJacobian();

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = glSystem_;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glSystem_;

  Teuchos::ParameterList &printParams = nlParamsPtr_->sublist( "Printing" );

  Teuchos::ParameterList& lsParams = nlParamsPtr_
                                    ->sublist ( "Direction" )
                                     .sublist ( "Newton" )
                                     .sublist ( "Linear Solver" );

  // Get initial solution
  Teuchos::RCP<Epetra_Vector> soln = glSystem_->getSolution();

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
    Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( printParams,
                                                          lsParams,
                                                          iReq,
                                                          iJac,
                                                          Analytic,
                                                          *soln ) );

  // Create the Group
  NOX::Epetra::Vector initialGuess ( soln, NOX::Epetra::Vector::CreateView );

  grpPtr_ = Teuchos::rcp ( new NOX::Epetra::Group ( printParams,
                                                    iReq,
                                                    initialGuess,
                                                    linSys ) );
}
// =============================================================================
// Create the solver
void
glNox::createSolver()
{
  if ( grpPtr_.is_null() )
      throw glException( "glNox::createSolver",
                         "Group not initialized" );

  if ( combo_.is_null() )
      throw glException( "glNox::createSolver",
                         "Combo not initialized" );

  if ( nlParamsPtr_.is_null() )
      throw glException( "glNox::createSolver",
                         "Nonlinear solver parameters not initialized" );

  solver_ = NOX::Solver::buildSolver ( grpPtr_,
                                       combo_,
                                       nlParamsPtr_ );
}
// =============================================================================
// TODO:
// Look into having this done by Trilinos. If executed on a multiproc
// environment, we don't want p to be fully present on all processors.
void
glNox::reOrder( Tpetra::Vector<double_complex>    &psi,
                const Teuchos::RCP<StaggeredGrid> &sGrid )
{
  int NumElements = psi.getGlobalLength();

  // fill p:
  std::vector<int> p ( NumElements );

  // copy over
  Tpetra::Vector<double_complex,int> psiTmp( psi );

//   = Teuchos::rcp ( new Tpetra::Vector<double_complex,int> ( psiLexicographic->getMap(),1 ) );
  sGrid->lexicographic2grid ( &p );

  Teuchos::ArrayRCP<const double_complex> psiTmpView = psiTmp.get1dView();
  for ( int k=0; k<NumElements; k++ )
    {
      psi.replaceGlobalValue ( p[k],
	                       psiTmpView[k]
			     );
    }
}
// =============================================================================
void
glNox::setNonlinearSolverParameters( Teuchos::ParameterList & nlParams )
{

  // Set the nonlinear solver method
  nlParams.set ( "Nonlinear Solver", "Line Search Based" );

  Teuchos::ParameterList& printParams = nlParams.sublist ( "Printing" );
  setPrintParameters( printParams );

  Teuchos::ParameterList& searchParams = nlParams.sublist ( "Line Search" );
  setSearchParameters( searchParams );

  Teuchos::ParameterList& dirParams = nlParams.sublist ( "Direction" );
  setDirectionParameters( dirParams );

  // Let's force all status tests to do a full check
  nlParams.sublist ( "Solver Options" )
                                   .set( "Status Test Check Type", "Complete" );

  if ( verbose_ ) // get custom pre/post actions
    {
      Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo =
        Teuchos::rcp ( new GlPrePostOperator ( glSystem_,
                                               problemParameters_ ) );
      nlParamsPtr_->sublist ( "Solver Options" )
                                 .set ( "User Defined Pre/Post Operator", ppo );
    }
}
// =============================================================================
  // Set the printing parameters in the "Printing" sublist
void
glNox::setPrintParameters( Teuchos::ParameterList & printParams )
{
  printParams.set ( "MyPID", MyPID_ );
  printParams.set ( "Output Precision", 16 );
  printParams.set ( "Output Processor", 0 );

  if ( verbose_ )
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
}
// =============================================================================
void
glNox::setSearchParameters( Teuchos::ParameterList & searchParams )
{
  searchParams.set ( "Method", "Full Step" );
}
// =============================================================================
void
glNox::setDirectionParameters( Teuchos::ParameterList & dirParams )
{
  dirParams.set ( "Method", "Newton" );

  Teuchos::ParameterList& newtonParams = dirParams.sublist ( "Newton" );
  setNewtonParameters( newtonParams );
}
// =============================================================================
void
glNox::setNewtonParameters( Teuchos::ParameterList & newtonParams )
{
  newtonParams.set ( "Forcing Term Method", "Constant" );

  Teuchos::ParameterList& lsParams = newtonParams.sublist ( "Linear Solver" );
  setLinearSolverParameters( lsParams );
}
// =============================================================================
void
glNox::setLinearSolverParameters( Teuchos::ParameterList & lsParams )
{
//  lsParams.set("Amesos Solver", "Amesos_Superlu");

  // lsParams.set("Aztec Solver", "BiCGStab");
  lsParams.set ( "Aztec Solver", "GMRES" );
  lsParams.set ( "Output Frequency", 32 );
  lsParams.set ( "Output Solver Details", true );
  lsParams.set ( "Max Iterations", 2000 );
  lsParams.set ( "Tolerance", 1e-4 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
// =============================================================================
void
glNox::createConvergenceTests()
{
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
                        Teuchos::rcp ( new NOX::StatusTest::NormF ( 1.0e-10 ) );

  if ( grpPtr_.is_null() )
      throw glException( "glNox::createConvergenceTests",
                         "Group not initialized" );

  NOX::Epetra::Group& grp = *grpPtr_;

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

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp ( new NOX::StatusTest::MaxIters ( maxNonlinearIterations_ ) );
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp ( new NOX::StatusTest::FiniteValue );
  // this test is useful if we start in a solution
  Teuchos::RCP<NOX::StatusTest::NormF> absresexact =
    Teuchos::rcp ( new NOX::StatusTest::NormF ( 1.0e-13 ) );

  combo_ =
     Teuchos::rcp ( new NOX::StatusTest::Combo ( NOX::StatusTest::Combo::OR ) );

  combo_->addStatusTest ( fv );
  combo_->addStatusTest ( absresexact );
  combo_->addStatusTest ( converged );
  combo_->addStatusTest ( maxiters );
}
// =============================================================================
double
glNox::computeJacobianConditionNumber()
{
  double kappa;
  try
    {
      const NOX::Epetra::Group& finalGroup =
        dynamic_cast<const NOX::Epetra::Group&> ( solver_->getSolutionGroup() );
      grpPtr_->computeJacobian();
      grpPtr_->computeJacobianConditionNumber ( 2000, 1e-2, 30, true );
      kappa = finalGroup.getJacobianConditionNumber();
    }
  catch ( std::exception& e )
    {
      std::cerr << e.what() << std::endl;
    }

  return kappa;
}
// =============================================================================
// compute some eigenvalues using anasazi
void
glNox::computeJacobianEigenvalues()
{
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

  const NOX::Epetra::Group& finalGroup =
        dynamic_cast<const NOX::Epetra::Group&> ( solver_->getSolutionGroup() );

  const Epetra_Vector& finalSolution =
            ( dynamic_cast<const NOX::Epetra::Vector&> ( finalGroup.getX() ) ).
                                                              getEpetraVector();

  Epetra_BlockMap Map = finalSolution.Map();

  // shortname for the final Jacobian
  Teuchos::RCP<Epetra_Operator> J = 
                              grpPtr_->getLinearSystem()->getJacobianOperator();

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
       std::string message( "Anasazi::BasicEigenproblem::setProblem() returned with error." );
       throw glException( "glNox::computeJacobianEigenvalues", message );
    }

  // Initialize the Block Arnoldi solver
//       Anasazi::BlockDavidsonSolMgr<double, MV, OP> 
//       Anasazi::LOBPCGSolMgr<double, MV, OP> 
//       Anasazi::RTRSolMgr<double, MV, OP> 
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> 
					        MySolverMgr ( MyProblem, MyPL );

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if ( returnCode != Anasazi::Converged && MyPID_==0 && verbose_ )
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
// =============================================================================
// void
// glNox::getSolution()
// {
//   const NOX::Epetra::Group& finalGroup =
//         dynamic_cast<const NOX::Epetra::Group&> ( solver_->getSolutionGroup() );
// 
//   const Epetra_Vector& finalSolution =
//             ( dynamic_cast<const NOX::Epetra::Vector&> ( finalGroup.getX() ) ).
//                                                               getEpetraVector();
// }
// =============================================================================
  // print the solution to a file
void
glNox::printSolutionToFile( std::string fileName )
{
  const NOX::Epetra::Group& finalGroup =
	dynamic_cast<const NOX::Epetra::Group&> ( solver_->getSolutionGroup() );

  const Epetra_Vector& finalSolution =
	    ( dynamic_cast<const NOX::Epetra::Vector&> ( finalGroup.getX() ) ).
							      getEpetraVector();

  glSystem_->solutionToFile ( finalSolution,
                              problemParameters_,
                              fileName );
}
// =============================================================================
int
glNox::checkConvergence()
{
  int status = 0;

//   // 1. Convergence
//   if ( solvStatus != NOX::StatusTest::Converged ) {
//       status = 1;
//       if ( printing.isPrintType ( NOX::Utils::Error ) )
//         printing.out() << "Nonlinear solver failed to converge!" << endl;
//   }

#ifndef HAVE_MPI
  // 2. Linear solve iterations (53) - SERIAL TEST ONLY!
  //    The number of linear iterations changes with # of procs.
  if ( const_cast<Teuchos::ParameterList&> ( solver_->getList() )
                                           .sublist ( "Direction" )
                                           .sublist ( "Newton" )
                                           .sublist ( "Linear Solver" )
                                           .sublist ( "Output" )
                                           .get ( "Total Number of Linear Iterations",0 )
     != 53 )
    {
      status = 2;
    }
#endif
  // 3. Nonlinear solve iterations (10)
  if ( const_cast<Teuchos::ParameterList&> ( solver_->getList() )
                                             .sublist ( "Output" )
                                             .get ( "Nonlinear Iterations", 0 )
        == maxNonlinearIterations_ )
    status = 3;

  return status;
}
// =============================================================================
void
glNox::setVerbose( bool verbose )
{
    verbose_ = verbose;
}
// =============================================================================