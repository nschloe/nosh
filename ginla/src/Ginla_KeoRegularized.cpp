// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2010, 2011  Nico Schl\"omer
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
typedef double                           ST;
typedef Epetra_MultiVector               MV;
typedef Epetra_Operator                  OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
// =============================================================================
namespace Ginla {
// =============================================================================
KeoRegularized::
KeoRegularized( const Teuchos::RCP<Ginla::KeoFactory> & keoFactory
                 ):
        useTranspose_ ( false ),
        comm_( Teuchos::rcpFromRef( keoFactory->getMesh()->getComm() ) ),
        keoFactory_( keoFactory ),
        keoRegularized_( Teuchos::null ),
        keoMlPrec_ ( Teuchos::null ),
        keoIluProblem_( Teuchos::null ),
        keoIluSolver_( Teuchos::null ),
        invType_( INVERT_ML ),
#ifdef GINLA_TEUCHOS_TIME_MONITOR
        timerRebuild_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoRegularized::rebuild") ),
        timerRebuildMl_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoRegularized::rebuildMl") ),
        timerRebuildIlu_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoRegularized::rebuildIlu") ),
        timerRegularization_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoRegularized::rebuild:regularization") ),
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
SetUseTranspose( bool UseTranspose )
{
    useTranspose_ = UseTranspose;
    return 0;
}
// =============================================================================
int
KeoRegularized::
Apply ( const Epetra_MultiVector & X,
              Epetra_MultiVector & Y
      ) const
{
    TEUCHOS_ASSERT( !keoRegularized_.is_null() );
    return keoRegularized_->Apply( X, Y );
}
// =============================================================================
int
KeoRegularized::
ApplyInverse( const Epetra_MultiVector & X,
                    Epetra_MultiVector & Y
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
        TEST_FOR_EXCEPT_MSG( true,
                                     "Illegal value of the invType \""
                                     << invType_ << "\"." );
        break;
  }

  return err;
}
// =============================================================================
int
KeoRegularized::
ApplyInverseMl_( const Epetra_MultiVector & X,
                       Epetra_MultiVector & Y
               ) const
{
   bool verbose = false;
   int frequency = 10;

   // Belos part
   Teuchos::ParameterList belosList;
   belosList.set( "Convergence Tolerance", 1.0e-10 ); // Relative convergence tolerance requested
   if (verbose) {
     belosList.set( "Verbosity",
                    Belos::Errors +
                    Belos::Warnings +
                    Belos::TimingDetails +
                    Belos::StatusTestDetails
                  );
     if (frequency > 0)
       belosList.set( "Output Frequency", frequency );
   }
   else
     belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

   // Make sure to have a solid initial guess.
   // Belos, for example, does not initialze Y before passing it here.
   Y.PutScalar( 0.0 );

   TEUCHOS_ASSERT( !keoRegularized_.is_null() );

   // Construct an unpreconditioned linear problem instance.
   Teuchos::RCP<const Epetra_MultiVector> Xptr = Teuchos::rcpFromRef( X );
   Teuchos::RCP<Epetra_MultiVector> Yptr = Teuchos::rcpFromRef( Y );
   Belos::LinearProblem<double,MV,OP> problem( keoRegularized_, Yptr, Xptr );
   // Make sure the problem sets up correctly.
   TEUCHOS_ASSERT( problem.setProblem() );
   // -------------------------------------------------------------------------
   // add preconditioner
   TEUCHOS_ASSERT( !keoMlPrec_.is_null() );
   problem.setLeftPrec( keoMlPrec_ );
   // -------------------------------------------------------------------------
   // Create an iterative solver manager.
   Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver =
           Teuchos::rcp( new Belos::PseudoBlockCGSolMgr<double,MV,OP>( Teuchos::rcp(&problem,false),
                                                                       Teuchos::rcp(&belosList,false)
                                                                     )
                       );

   // Perform solve
   Belos::ReturnType ret = newSolver->solve();

   /* // compute the residual explicitly
   Epetra_MultiVector R( Y.Map(), Y.NumVectors() );
   keoRegularized_->Apply( Y, R );
   R.Update( 1.0, X, -1.0 );
   double s[1];
   R.Norm2( s );
   std::cout << " PRECON RES = " << s[0] << " in " << newSolver->getNumIters() << " iters" << std::endl;*/

   return ret==Belos::Converged ? 0 : -1;
}
// =============================================================================
int
KeoRegularized::
ApplyInverseIlu_ ( const Epetra_MultiVector & X,
                         Epetra_MultiVector & Y
                 ) const
{
    TEUCHOS_ASSERT( !keoIluProblem_.is_null() );
    TEUCHOS_ASSERT( !keoIluSolver_.is_null() );

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
NormInf () const
{
    TEST_FOR_EXCEPT( "Not yet implemented." );
    return 0.0;
}
// =============================================================================
const char *
KeoRegularized::
Label () const
{
    return "Kinetic energy operator for Ginzburg--Landau";
}
// =============================================================================
bool
KeoRegularized::
UseTranspose () const
{
    return useTranspose_;
}
// =============================================================================
bool
KeoRegularized::
HasNormInf () const
{
    return false;
}
// =============================================================================
const Epetra_Comm &
KeoRegularized::
Comm () const
{
    TEUCHOS_ASSERT( !comm_.is_null() );
    return *comm_;
}
// =============================================================================
const Epetra_Map &
KeoRegularized::
OperatorDomainMap () const
{
    TEUCHOS_ASSERT( !keoRegularized_.is_null() );
    return keoRegularized_->OperatorDomainMap();
}
// =============================================================================
const Epetra_Map &
KeoRegularized::
OperatorRangeMap () const
{
    TEUCHOS_ASSERT( !keoRegularized_.is_null() );
    return keoRegularized_->OperatorRangeMap();
}
// =============================================================================
void
KeoRegularized::
rebuild()
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm(*timerRebuild_);
#endif
    // -------------------------------------------------------------------------
    // rebuild the keo
    keoRegularized_ = keoFactory_->buildKeo();
    keoRegularized_->Scale( -1.0 );
    // -------------------------------------------------------------------------
    // regularization
    double mu = keoFactory_->getMvpParameters()->getValue( "mu" );
    if ( fabs( mu ) < 1.0e-5 )
    {
#ifdef GINLA_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor tm(*timerRegularization_);
#endif
        // Add a regularization to the diagonal.
        Epetra_Vector e( keoRegularized_->DomainMap() );
        TEUCHOS_ASSERT_EQUALITY( 0, e.PutScalar( 1.0 ) );

        Epetra_Vector diag( keoRegularized_->DomainMap() );
        TEUCHOS_ASSERT_EQUALITY( 0, keoRegularized_->ExtractDiagonalCopy( diag ) );
        TEUCHOS_ASSERT_EQUALITY( 0, diag.Update( 1.0e-3, e, 1.0 ) );
        TEUCHOS_ASSERT_EQUALITY( 0, keoRegularized_->ReplaceDiagonalValues( diag ) );
    }
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
          TEST_FOR_EXCEPT_MSG( true,
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
rebuild( const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams
       )
{
    keoFactory_->updateParameters( mvpParams );
    this->rebuild();
    return;
}
// =============================================================================
void
KeoRegularized::
rebuildMl_()
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm(*timerRebuildMl_);
#endif
    // rebuild ML structure
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults( "SA", MLList );
//     MLList.set("ML output", 0);
    MLList.set("max levels", 10);
    MLList.set("increasing or decreasing", "increasing");
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("smoother: type", "Chebyshev"); // "block Gauss-Seidel" "Chebyshev"
//     MLList.set("aggregation: threshold", 0.0);
    MLList.set("smoother: sweeps", 3);
    MLList.set("smoother: pre or post", "both");
    MLList.set("coarse: type", "Amesos-KLU");
    MLList.set("PDE equations", 2);

    TEUCHOS_ASSERT( !keoRegularized_.is_null() );

    Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec =
        Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(*keoRegularized_, MLList) );
    TEUCHOS_ASSERT( !MLPrec.is_null() );
//     MLPrec->PrintUnused(0);

    // Create the Belos preconditioned operator from the preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    keoMlPrec_ = Teuchos::rcp( new Belos::EpetraPrecOp( MLPrec ) );

    return;
}
// =============================================================================
void
KeoRegularized::
rebuildIlu_()
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm(*timerRebuildIlu_);
#endif
    // set the matrix the linear problem
    TEUCHOS_ASSERT( !keoRegularized_.is_null() );
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