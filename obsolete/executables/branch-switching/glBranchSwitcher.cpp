#include "glBranchSwitcher.h"

#include "glException.h"

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

// =============================================================================                       
GlBranchSwitcher::GlBranchSwitcher ( const Teuchos::RCP<GlSystem> glSystem,
	                             const double ds,
                                     const double alpha0,
                                     const double alpham1,
                                     const Teuchos::RCP<const Epetra_Vector> u0,
                                     const Teuchos::RCP<const Epetra_Vector> um1 ) :
  glSystem_(glSystem),
  ds_(ds),
  alpha0_(alpha0),
  alpham1_(alpham1),
  u0_( u0 ),
  um1_( um1 ),
  dFdalpha_( 0 )
{

  computePhi2();

}
// =============================================================================                       
GlBranchSwitcher::~GlBranchSwitcher()
{
}
// =============================================================================                       
bool
GlBranchSwitcher::computeF ( const Epetra_Vector &x,
                             Epetra_Vector       &F,
                             const NOX::Epetra::Interface::Required::FillType fillFlag
         )
{
}
// =============================================================================                       
bool
GlBranchSwitcher::computeJacobian ( const Epetra_Vector &x,
                  Epetra_Operator     &Jac )
{
}
// =============================================================================                       
bool
GlBranchSwitcher::computePreconditioner ( const Epetra_Vector     &x,
                                          Epetra_Operator         &Prec,
                                          Teuchos::ParameterList* precParams
                                        )  const
{
}
// =============================================================================                       
Teuchos::RCP<Epetra_Vector>
GlBranchSwitcher::getSolution() const
{
}
// =============================================================================                       
Teuchos::RCP<Epetra_CrsMatrix>
GlBranchSwitcher::getJacobian() const
{
}
// =============================================================================                       
void
GlBranchSwitcher::computePhi2()
{
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
GlBranchSwitcher::getNullvector( const Teuchos::RCP<const Epetra_Operator> A )
{
  bool verbose = true;
  bool debug = true;
  std::string which ( "SM" ); // smallest magnitude

  bool boolret;

  typedef double ScalarType;
  typedef Teuchos::ScalarTraits<ScalarType>          SCT;
  typedef SCT::magnitudeType                         MagnitudeType;
  typedef Epetra_MultiVector                         MV;
  typedef Epetra_Operator                            OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;

  Epetra_BlockMap Map = u0_->Map();

  // - - - - - - - - - - - - - - - - -
  // Start the block Arnoldi iteration
  // - - - - - - - - - - - - - - - - -
  // Variables used for the Block Krylov Schur Method
  int nev = 1; // only get *one* eigenvector here (the one corresponding
               // to the 0 eigenvalue
  int blockSize = 10;
  int numBlocks = 10;
  int maxRestarts = 25;
  double tol = 1.0e-10;
  
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
    Teuchos::rcp ( new Anasazi::BasicEigenproblem<double, MV, OP> ( A, ivec ) );

  // Set the number of eigenvalues requested
  MyProblem->setNEV ( nev );

  int MyPID = Map.Comm().MyPID();

  // Inform the eigenproblem that you are finishing passing it information
  bool ierr = MyProblem->setProblem();
  if ( ierr != true )
  {
  if ( MyPID == 0 )
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
#ifdef HAVE_MPI
      MPI_Finalize() ;
#endif
      return Teuchos::null;
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

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<ScalarType> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  std::vector<int> index = sol.index;
  int numev = sol.numVecs;

  double resTol = 1.0e-14;
  double realPart = evals[0].realpart;
  double imagPart = evals[0].imagpart;
  double abs0 = sqrt( realPart*realPart + imagPart*imagPart );
  if ( abs0>resTol ) {
      std::string message = "Computes eigenvalue not 0";
      throw glException( "GlBranchSwitcher::computeNullvector", message );
  }

  Teuchos::RCP<Epetra_Vector> ev1 ( ( *evecs ) ( 0 ) );

  return ev1;

}
// =============================================================================
void
GlBranchSwitcher::computedFdalpha()
{

  LOCA::ParameterVector paramVec;
  int ierr;

  // set parameter
  Teuchos::RCP<Epetra_Vector>  FAlpha0 = Teuchos::rcp( new Epetra_Vector(*u0_) );
  paramVec.setValue ( "H0", alpha0_ );
  glSystem_->setParameters ( paramVec );
  ierr = glSystem_->computeF ( *u0_, *FAlpha0 );

  Teuchos::RCP<Epetra_Vector>  dFdalpha = Teuchos::rcp( new Epetra_Vector(*u0_) );
  paramVec.setValue ( "H0", alpham1_ );
  glSystem_->setParameters ( paramVec );
  ierr = glSystem_->computeF ( *um1_, *dFdalpha );

  // compute finite difference
  double dalpha = alpha0_-alpham1_;
  dFdalpha->Update( 1.0/dalpha, *FAlpha0, -1.0/dalpha );
  
}
// =============================================================================                       
