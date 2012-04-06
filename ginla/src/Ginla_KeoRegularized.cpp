// @HEADER
//
//    Regularized kinetic energy operator.
//    Copyright (C) 2010--2012  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER

#include "Ginla_KeoRegularized.hpp"
#include "Ginla_KeoFactory.hpp"
#include "Ginla_StkMesh.hpp"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <ml_epetra_preconditioner.h>

#include <Amesos_BaseSolver.h>
#include <Epetra_LinearProblem.h>

#include <Teuchos_VerboseObject.hpp>

#include <Amesos.h>

// =============================================================================
// some typdefs for Belos
typedef double ST;
typedef Epetra_MultiVector MV;
typedef Epetra_Operator OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
// =============================================================================
namespace Ginla {
// =============================================================================
KeoRegularized::
KeoRegularized( const Teuchos::RCP<const Ginla::StkMesh> &mesh,
                const Teuchos::RCP<const Epetra_Vector> &thickness,
                const Teuchos::RCP<Ginla::KeoFactory> &keoFactory,
                const Teuchos::RCP<const Epetra_Vector> &psi
              ) :
  useTranspose_( false ),
  mesh_( mesh ),
  thickness_( thickness ),
  comm_( keoFactory->getComm() ),
  keoFactory_( keoFactory ),
  psi_( psi ),
  keoRegularized_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, keoFactory_->getKeo()->Graph()))),
  MlPrec_( Teuchos::null ),
  keoIluProblem_( Teuchos::null ),
  keoIluSolver_( Teuchos::null ),
  invType_( INVERT_ML ),
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  timerRebuild0_( Teuchos::TimeMonitor::getNewTimer(
                   "Ginla: KeoRegularized::rebuild::ML init" ) ),
  timerRebuild1_( Teuchos::TimeMonitor::getNewTimer(
                   "Ginla: KeoRegularized::rebuild::ML rebuild" ) ),
#endif
  out_( Teuchos::VerboseObjectBase::getDefaultOStream() )
{
}
// =============================================================================
KeoRegularized::
~KeoRegularized()
{
}
// =============================================================================
int
KeoRegularized::
SetUseTranspose( bool useTranspose )
{
  useTranspose_ = useTranspose;
  return 0;
}
// =============================================================================
int
KeoRegularized::
Apply( const Epetra_MultiVector &X,
       Epetra_MultiVector &Y
       ) const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keoRegularized_.is_null() );
#endif
  return keoRegularized_->Apply( X, Y );
}
// =============================================================================
int
KeoRegularized::
ApplyInverse( const Epetra_MultiVector &X,
              Epetra_MultiVector &Y
              ) const
{
  int err;

  switch ( invType_ )
  {
    case INVERT_ILU:
      err = this->ApplyInverseIlu_( X, Y );
      break;
    case INVERT_ML:
      err = this->ApplyInverseMl_( X, Y );
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                           "Illegal value of the invType \""
                           << invType_ << "\"." );
      break;
  }

  return err;
}
// =============================================================================
int
KeoRegularized::
ApplyInverseMl_(const Epetra_MultiVector &X,
                Epetra_MultiVector &Y
                ) const
{
  // Just apply one (inverse) AMG cycle.
  return MlPrec_->ApplyInverse(X, Y);

//   bool verbose = false;
//   int frequency = 10;
// 
//   // Belos part
//   Teuchos::ParameterList belosList;
//   // Relative convergence tolerance requested
//   // TODO This could be replaced by something adaptive in the future.
//   belosList.set( "Convergence Tolerance", 1.0e-10 );
//   if (verbose) {
//     belosList.set( "Verbosity",
//                    Belos::Errors +
//                    Belos::Warnings +
//                    Belos::TimingDetails +
//                    Belos::StatusTestDetails
//                    );
//     if (frequency > 0)
//       belosList.set( "Output Frequency", frequency );
//   }
//   else
//     belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
// 
//   // Make sure to have a solid initial guess.
//   // Belos, for example, does not initialze Y before passing it here.
//   Y.PutScalar( 0.0 );
// 
// #ifdef _DEBUG_
//   TEUCHOS_ASSERT( !keoRegularized_.is_null() );
// #endif
// 
//   // Construct an unpreconditioned linear problem instance.
//   Teuchos::RCP<const Epetra_MultiVector> Xptr = Teuchos::rcpFromRef( X );
//   Teuchos::RCP<Epetra_MultiVector> Yptr = Teuchos::rcpFromRef( Y );
//   Belos::LinearProblem<double,MV,OP> problem( keoRegularized_, Yptr, Xptr );
//   // Make sure the problem sets up correctly.
//   TEUCHOS_ASSERT( problem.setProblem() );
//   // -------------------------------------------------------------------------
//   // add preconditioner
//   // Create the Belos preconditioned operator from the preconditioner.
//   // NOTE:  This is necessary because Belos expects an operator to apply the
//   //        preconditioner with Apply() NOT ApplyInverse().
//   Teuchos::RCP<Belos::EpetraPrecOp> keoMlPrec =
//       Teuchos::rcp( new Belos::EpetraPrecOp( MlPrec_ ) );
//   problem.setLeftPrec( keoMlPrec );
//   // -------------------------------------------------------------------------
//   // Create an iterative solver manager.
//   Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver =
//     Teuchos::rcp( new Belos::PseudoBlockCGSolMgr<double,MV,
//                                                  OP>( Teuchos::rcp( &problem,
//                                                                     false ),
//                                                       Teuchos::rcp( &
//                                                                     belosList,
//                                                                     false )
//                                                       )
//                   );
// 
//   // Perform solve
//   Belos::ReturnType ret = newSolver->solve();
// 
// //    // compute the residual explicitly
// //    Epetra_MultiVector R( Y.Map(), Y.NumVectors() );
// //    keoRegularized_->Apply( Y, R );
// //    R.Update( 1.0, X, -1.0 );
// //    double s[1];
// //    R.Norm2( s );
// //    std::cout << " PRECON RES = " << s[0] << " in " << newSolver->getNumIters() << " iters" << std::endl;
// 
//   return ret==Belos::Converged ? 0 : -1;
}
// =============================================================================
int
KeoRegularized::
ApplyInverseIlu_( const Epetra_MultiVector &X,
                  Epetra_MultiVector &Y
                  ) const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keoIluProblem_.is_null() );
  TEUCHOS_ASSERT( !keoIluSolver_.is_null() );
#endif

  // set left- and right-hand side
  keoIluProblem_->SetLHS( &Y );
  Epetra_MultiVector T = X;
  keoIluProblem_->SetRHS( &T );

  // solve and return error code
  return keoIluSolver_->Solve();
}
// =============================================================================
double
KeoRegularized::
NormInf() const
{
  TEUCHOS_TEST_FOR_EXCEPT( "Not yet implemented." );
  return 0.0;
}
// =============================================================================
const char *
KeoRegularized::
Label() const
{
  return "Kinetic energy operator for Ginzburg--Landau";
}
// =============================================================================
bool
KeoRegularized::
UseTranspose() const
{
  return useTranspose_;
}
// =============================================================================
bool
KeoRegularized::
HasNormInf() const
{
  return false;
}
// =============================================================================
const Epetra_Comm &
KeoRegularized::
Comm() const
{
  return comm_;
}
// =============================================================================
const Epetra_Map &
KeoRegularized::
OperatorDomainMap() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keoRegularized_.is_null() );
#endif
  return keoRegularized_->OperatorDomainMap();
}
// =============================================================================
const Epetra_Map &
KeoRegularized::
OperatorRangeMap() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keoRegularized_.is_null() );
#endif
  return keoRegularized_->OperatorRangeMap();
}
// =============================================================================
void
KeoRegularized::
rebuild(const Teuchos::RCP<const Epetra_Vector> & psi)
{
  psi_ = psi;
  // -------------------------------------------------------------------------
  // Copy over the matrix.
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keoFactory_.is_null() );
#endif
  *keoRegularized_ = *(keoFactory_->getKeo());
  keoRegularized_->Scale( -1.0 );
  // -------------------------------------------------------------------------
  // Add 2*|psi|^2 to the diagonal.
  Epetra_Vector diag( keoRegularized_->RowMap() );
  TEUCHOS_ASSERT_EQUALITY( 0, keoRegularized_->ExtractDiagonalCopy( diag ) );
  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());
  int numMyPoints = controlVolumes.MyLength();
  for ( int k=0; k<numMyPoints; k++ )
  {
    double alpha = 2.0 * controlVolumes[k] * (*thickness_)[k]
                 * ((*psi_)[2*k]*(*psi_)[2*k] + (*psi_)[2*k+1]*(*psi_)[2*k+1]);
    diag[2*k] += alpha;
    diag[2*k+1] += alpha;
  }
  TEUCHOS_ASSERT_EQUALITY( 0, keoRegularized_->ReplaceDiagonalValues( diag ) );
  // -------------------------------------------------------------------------
  // Rebuild preconditioner for this object. Not to be mistaken for the
  // object itself.
  switch ( invType_ )
  {
    case INVERT_ILU:
      return this->rebuildIlu_();
      break;
    case INVERT_ML:
      return this->rebuildMl_();
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                           "Illegal value of the invType \""
                           << invType_ << "\"."
                           );
      break;
  }
  // -------------------------------------------------------------------------

  return;
}
// =============================================================================
void
KeoRegularized::
rebuild(const Teuchos::RCP<const LOCA::ParameterVector> &mvpParams,
        const Teuchos::RCP<const Epetra_Vector> & psi
        )
{
  keoFactory_->updateParameters( mvpParams );
  this->rebuild( psi );
  return;
}
// =============================================================================
void
KeoRegularized::
rebuildMl_()
{
  // For reusing the ML structure, see
  // http://trilinos.sandia.gov/packages/docs/dev/packages/ml/doc/html/classML__Epetra_1_1MultiLevelPreconditioner.html#a0a5c1d47c6938d2ec1cb9bb710723c1e
  if ( MlPrec_.is_null() )
  {
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm( *timerRebuild0_ );
#endif
    // build ML structure
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults( "SA", MLList );
//     MLList.set("ML output", 0);
    MLList.set( "max levels", 10 );
    MLList.set( "increasing or decreasing", "increasing" );
    MLList.set( "aggregation: type", "Uncoupled" );
    MLList.set( "smoother: type", "Chebyshev" );     // "block Gauss-Seidel" "Chebyshev"
    MLList.set( "aggregation: threshold", 0.0 );
    MLList.set( "smoother: sweeps", 3 );
    MLList.set( "smoother: pre or post", "both" );
    MLList.set( "coarse: type", "Amesos-KLU" );
    MLList.set( "PDE equations", 2 );
    // reuse the multilevel hierarchy
    MLList.set("reuse: enable", true);

#ifdef _DEBUG_
    TEUCHOS_ASSERT( !keoRegularized_.is_null() );
#endif

    // From http://trilinos.sandia.gov/packages/docs/r10.8/packages/ml/doc/html/classML__Epetra_1_1MultiLevelPreconditioner.html:
    // "It is important to note that ML is more restrictive than Epetra for
    //  the definition of maps. It is required that RowMatrixRowMap() is
    //  equal to OperatorRangeMap(). This is because ML needs to perform
    //  matrix-std::vector product, as well as getrow() functions, on the
    //  same data distribution.
    // "Also, for square matrices, OperatorDomainMap() must be as
    //  OperatorRangeMap()."
    // Make sure this is indeed the case.
#ifdef _DEBUG_
    TEUCHOS_ASSERT( keoRegularized_->OperatorRangeMap().
                    SameAs( keoRegularized_->RowMatrixRowMap() )
                  );
    TEUCHOS_ASSERT( keoRegularized_->OperatorDomainMap().
                    SameAs( keoRegularized_->OperatorRangeMap() )
                  );
#endif
    MlPrec_ =
      Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner( *keoRegularized_,
                                                            MLList ) );
  }
  else
  {
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm( *timerRebuild1_ );
#endif

    bool checkFiltering = true;
    TEUCHOS_ASSERT_EQUALITY( 0, MlPrec_->ComputePreconditioner(checkFiltering) );
    //TEUCHOS_ASSERT_EQUALITY( 0, MlPrec_->ReComputePreconditioner() );
  }

//    MlPrec_->PrintUnused(0);

  return;
}
// =============================================================================
void
KeoRegularized::
rebuildIlu_()
{
  // set the matrix the linear problem
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keoRegularized_.is_null() );
#endif
  if ( keoIluProblem_.is_null() )
    keoIluProblem_ = Teuchos::rcp( new Epetra_LinearProblem() );
  keoIluProblem_->SetOperator( &*keoRegularized_ );

  // do the factorization
  Amesos Factory;
  keoIluSolver_ = Teuchos::rcp( Factory.Create( "Klu", *keoIluProblem_ ) );

  // do symbolic and numerical factorizations
  // TODO reuse symbolic factorization
  TEUCHOS_ASSERT_EQUALITY( 0, keoIluSolver_->SymbolicFactorization() );
  TEUCHOS_ASSERT_EQUALITY( 0, keoIluSolver_->NumericFactorization() );

  return;
}
// =============================================================================
} // namespace Ginla
