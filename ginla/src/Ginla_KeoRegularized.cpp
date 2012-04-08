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
#include "Ginla_KeoContainer.hpp"
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
                const Teuchos::RCP<Ginla::KeoContainer> &keoContainer,
                const Teuchos::RCP<const Epetra_Vector> &psi
              ) :
  useTranspose_( false ),
  mesh_( mesh ),
  thickness_( thickness ),
  keoContainer_( keoContainer ),
  absPsiSquared_(Teuchos::rcp(new Epetra_Vector(*mesh->getComplexNonOverlapMap()))),
  keoRegularizedMatrix_( *keoContainer_->getKeo() ), // initialize the matrix with something (it gets copied over later anyways)
  comm_( keoContainer->getComm() ),
  MlPrec_( Teuchos::null ),
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  timerRebuild0_( Teuchos::TimeMonitor::getNewTimer(
                  "Ginla: KeoRegularized::rebuild::ML init" ) ),
  timerRebuild1_( Teuchos::TimeMonitor::getNewTimer(
                  "Ginla: KeoRegularized::rebuild::ML rebuild" ) ),
#endif
  out_( Teuchos::VerboseObjectBase::getDefaultOStream() )
{
  this->updateAbsPsiSquared_(psi);
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
  TEUCHOS_ASSERT( keoContainer_->getKeo()->DomainMap().SameAs( X.Map() ) );
  TEUCHOS_ASSERT( keoContainer_->getKeo()->RangeMap().SameAs( Y.Map() ) );
  TEUCHOS_ASSERT( Y.Map().SameAs( X.Map() ) );
  TEUCHOS_ASSERT( Y.Map().SameAs(absPsiSquared_->Map()) );
#endif
  // K * X ...
  TEUCHOS_ASSERT_EQUALITY(0, keoContainer_->getKeo()->Apply(X, Y));

  // - (K*X) + 2*|psi|*X.
  TEUCHOS_ASSERT_EQUALITY(0, Y.Multiply(2.0, *absPsiSquared_, X, -1.0));

  return 0;
}
// =============================================================================
int
KeoRegularized::
ApplyInverse(const Epetra_MultiVector &X,
             Epetra_MultiVector &Y
             ) const
{
  // Just apply one (inverse) AMG cycle.
  return MlPrec_->ApplyInverse(X, Y);

//   bool verbose = false;
// 
//   // Belos part
//   Teuchos::ParameterList belosList;
//   // Relative convergence tolerance requested.
//   // Set this to 0 and adapt the maximum number of iterations. This way,
//   // the preconditioner always does exactly the same thing (namely maxIter
//   // PCG iterations) and is independent of X. This avoids mathematical
//   // difficulties.
//   belosList.set( "Convergence Tolerance", 0.0 );
//   belosList.set( "Maximum Iterations", 1 );
//   if (verbose) {
//     belosList.set( "Verbosity",
//                    Belos::Errors +
//                    Belos::Warnings +
//                    Belos::TimingDetails +
//                    Belos::StatusTestDetails
//                    );
//     belosList.set( "Output Frequency", 10 );
//   }
//   else
//     belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
// 
//   // Make sure to have a solid initial guess.
//   // Belos, for example, does not initialze Y before passing it here.
//   Y.PutScalar( 0.0 );
// 
//   // Construct an unpreconditioned linear problem instance.
//   Teuchos::RCP<const Epetra_MultiVector> Xptr = Teuchos::rcpFromRef( X );
//   Teuchos::RCP<Epetra_MultiVector> Yptr = Teuchos::rcpFromRef( Y );
//   Teuchos::RCP<const Epetra_CrsMatrix> Aptr = Teuchos::rcpFromRef(keoRegularizedMatrix_);
//   Belos::LinearProblem<double,MV,OP> problem(Aptr, Yptr, Xptr);
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
//     Teuchos::rcp( new Belos::PseudoBlockCGSolMgr<double,MV,OP>
//       (Teuchos::rcp(&problem, false), Teuchos::rcp(&belosList, false))
//     );
// 
//   // Perform solve
//   Belos::ReturnType ret = newSolver->solve();
// 
// //    // compute the residual explicitly
// //    Epetra_MultiVector R( Y.Map(), Y.NumVectors() );
// //    keoRegularizedMatrix_.Apply( Y, R );
// //    R.Update( 1.0, X, -1.0 );
// //    double s[1];
// //    R.Norm2( s );
// //    std::cout << " PRECON RES = " << s[0] << " in " << newSolver->getNumIters() << " iters" << std::endl;
// 
// //   return ret==Belos::Converged ? 0 : -1;
//   return 0;
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
  return keoContainer_->getKeo()->OperatorDomainMap();
}
// =============================================================================
const Epetra_Map &
KeoRegularized::
OperatorRangeMap() const
{
  return keoContainer_->getKeo()->OperatorRangeMap();
}
// =============================================================================
void
KeoRegularized::
rebuildInverse()
{
  // -------------------------------------------------------------------------
  // Copy over the matrix.
  // This is necessary as we don't apply AMG to K, but to K+2|psi|^2.
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keoContainer_.is_null() );
  TEUCHOS_ASSERT( !keoContainer_->getKeo().is_null() );
#endif
  keoRegularizedMatrix_ = *keoContainer_->getKeo();
  keoRegularizedMatrix_.Scale( -1.0 );
  // -------------------------------------------------------------------------
  // Add 2*|psi|^2 to the diagonal.
  Epetra_Vector diag(keoRegularizedMatrix_.RowMap());
  TEUCHOS_ASSERT_EQUALITY(0, keoRegularizedMatrix_.ExtractDiagonalCopy(diag));
  diag.Update(2.0, *absPsiSquared_, 1.0);
  TEUCHOS_ASSERT_EQUALITY(0, keoRegularizedMatrix_.ReplaceDiagonalValues(diag));
  // -------------------------------------------------------------------------
  // Rebuild preconditioner for this object. Not to be mistaken for the
  // object itself.
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
    TEUCHOS_ASSERT( keoRegularizedMatrix_.OperatorRangeMap().
                    SameAs( keoRegularizedMatrix_.RowMatrixRowMap() )
                  );
    TEUCHOS_ASSERT( keoRegularizedMatrix_.OperatorDomainMap().
                    SameAs( keoRegularizedMatrix_.OperatorRangeMap() )
                  );
#endif
    MlPrec_ =
      Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(keoRegularizedMatrix_,
                                                           MLList));
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
  // -------------------------------------------------------------------------

  return;
}
// =============================================================================
void
KeoRegularized::
rebuildInverse(const Teuchos::RCP<const LOCA::ParameterVector> &mvpParams,
               const Teuchos::RCP<const Epetra_Vector> & psi
              )
{
  keoContainer_->updateParameters( mvpParams );
  this->updateAbsPsiSquared_(psi);
  this->rebuildInverse();
  return;
}
// =============================================================================
void
KeoRegularized::
updateAbsPsiSquared_(const Teuchos::RCP<const Epetra_Vector> &psi)
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !thickness_.is_null() );
  TEUCHOS_ASSERT( !psi.is_null() );
  TEUCHOS_ASSERT( !absPsiSquared_.is_null() );
#endif

  const Teuchos::RCP<const Epetra_Vector> &controlVolumes =
      mesh_->getControlVolumes();
  int numMyPoints = controlVolumes->MyLength();
  for ( int k=0; k<numMyPoints; k++ )
  {
    double alpha = (*controlVolumes)[k] * (*thickness_)[k]
                 * ((*psi)[2*k]*(*psi)[2*k] + (*psi)[2*k+1]*(*psi)[2*k+1]);
    (*absPsiSquared_)[2*k] = alpha;
    (*absPsiSquared_)[2*k+1] = alpha;
  }
  return;
}
// =============================================================================
} // namespace Ginla
