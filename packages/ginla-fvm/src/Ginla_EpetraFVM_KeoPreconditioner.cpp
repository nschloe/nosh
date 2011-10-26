/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "Ginla_EpetraFVM_KeoPreconditioner.hpp"
#include "Ginla_EpetraFVM_KeoFactory.hpp"
#include "Ginla_EpetraFVM_StkMesh.hpp"

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

// =============================================================================
// some typdefs for Belos
typedef double                           ST;
typedef Epetra_MultiVector               MV;
typedef Epetra_Operator                  OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
// =============================================================================
Ginla::EpetraFVM::KeoPreconditioner::
KeoPreconditioner( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh>               & mesh,
                   const Teuchos::RCP<const Epetra_Vector>                     & thickness,
                   const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp
                 ):
        useTranspose_ ( false ),
        comm_( Teuchos::rcpFromRef(mesh->getComm() ) ),
        keoFactory_( Teuchos::rcp( new Ginla::EpetraFVM::KeoFactory(mesh, thickness, mvp) ) ),
        keoPrec_( Teuchos::rcp( new Epetra_FECrsMatrix( Copy, keoFactory_->buildKeoGraph() ) ) ),
        belosPrec_ ( Teuchos::null ),
        keoProblem_( Teuchos::rcp( new Epetra_LinearProblem() ) ),
        keoSolver_( Teuchos::null ),
        invType_( Ml ),
        out_( Teuchos::VerboseObjectBase::getDefaultOStream() ),
        rebuildTime_( Teuchos::TimeMonitor::getNewTimer("KeoPreconditioner::rebuild") ),
        rebuildMlTime_( Teuchos::TimeMonitor::getNewTimer("KeoPreconditioner::rebuildMl") ),
        rebuildIluTime_( Teuchos::TimeMonitor::getNewTimer("KeoPreconditioner::rebuildIlu") ),
        regularizationTime_( Teuchos::TimeMonitor::getNewTimer("KeoPreconditioner::rebuild:regularization") )
{
}
// =============================================================================
Ginla::EpetraFVM::KeoPreconditioner::
~KeoPreconditioner()
{
}
// =============================================================================
int
Ginla::EpetraFVM::KeoPreconditioner::
SetUseTranspose( bool UseTranspose )
{
    useTranspose_ = UseTranspose;
    return 0;
}
// =============================================================================
int
Ginla::EpetraFVM::KeoPreconditioner::
Apply ( const Epetra_MultiVector & X,
              Epetra_MultiVector & Y
      ) const
{
    return keoPrec_->Apply( X, Y );
}
// =============================================================================
int
Ginla::EpetraFVM::KeoPreconditioner::
ApplyInverse( const Epetra_MultiVector & X,
                    Epetra_MultiVector & Y
            ) const
{
  int err;

  switch ( invType_ )
  {
    case Ilu:
        err = this->ApplyInverseIlu_( X, Y );
        break;
    case Ml:
        err = this->ApplyInverseMl_( X, Y );
        break;
    default:
        TEUCHOS_TEST_FOR_EXCEPTION( true,
                                    std::invalid_argument,
                                    "Illegal value of the invType \"" << invType_ << "\"."
                                  );
        break;
  }

  return err;
}
// =============================================================================
int
Ginla::EpetraFVM::KeoPreconditioner::
ApplyInverseMl_( const Epetra_MultiVector & X,
                       Epetra_MultiVector & Y
               ) const
{
   bool verbose = false;
   int frequency = 10;
   bool proc_verbose = verbose && (X.Comm().MyPID()==0);

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

   // check for NaNs in the initial guess
   double r[ Y.NumVectors() ];
   Y.Norm1( r );
   for ( int k=0; k<Y.NumVectors(); k++ )
       TEUCHOS_TEST_FOR_EXCEPTION( r[k]!=r[k],
                                   std::runtime_error,
                                   "The input guess appears to contain NaNs." );

   // Construct an unpreconditioned linear problem instance.
   Teuchos::RCP<const Epetra_MultiVector> Xptr = Teuchos::rcpFromRef( X );
   Teuchos::RCP<Epetra_MultiVector> Yptr = Teuchos::rcpFromRef( Y );
   Belos::LinearProblem<double,MV,OP> problem( keoPrec_, Yptr, Xptr );
   bool set = problem.setProblem();

   if ( !set )
   {
     *out_ << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
     return -1;
   }
   // -------------------------------------------------------------------------
   // add preconditioner
   TEUCHOS_ASSERT( !belosPrec_.is_null() );
   problem.setLeftPrec( belosPrec_ );
   // -------------------------------------------------------------------------
   // Create an iterative solver manager.
   Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver =
           Teuchos::rcp( new Belos::PseudoBlockCGSolMgr<double,MV,OP>( Teuchos::rcp(&problem,false),
                                                                       Teuchos::rcp(&belosList,false)
                                                                     )
                       );

   // Perform solve
   Belos::ReturnType ret = newSolver->solve();

   return ret==Belos::Converged ? 0 : -1;
}
// =============================================================================
int
Ginla::EpetraFVM::KeoPreconditioner::
ApplyInverseIlu_ ( const Epetra_MultiVector & X,
                         Epetra_MultiVector & Y
                 ) const
{
    TEUCHOS_ASSERT( !keoProblem_.is_null() );
    TEUCHOS_ASSERT( !keoSolver_.is_null() );

    // set left- and right-hand side
    keoProblem_->SetLHS( &Y );
    Epetra_MultiVector T = X;
    keoProblem_->SetRHS( &T );

    // solve and return error code
    return keoSolver_->Solve();
}
// =============================================================================
double
Ginla::EpetraFVM::KeoPreconditioner::
NormInf () const
{
    TEUCHOS_TEST_FOR_EXCEPTION( true,
                                std::logic_error,
                                "Not yet implemented." );
    return 0.0;
}
// =============================================================================
const char *
Ginla::EpetraFVM::KeoPreconditioner::
Label () const
{
    return "Kinetic energy operator for Ginzburg--Landau";
}
// =============================================================================
bool
Ginla::EpetraFVM::KeoPreconditioner::
UseTranspose () const
{
    return useTranspose_;
}
// =============================================================================
bool
Ginla::EpetraFVM::KeoPreconditioner::
HasNormInf () const
{
    return false;
}
// =============================================================================
const Epetra_Comm &
Ginla::EpetraFVM::KeoPreconditioner::
Comm () const
{
    TEUCHOS_ASSERT( !comm_.is_null() );
    return *comm_;
}
// =============================================================================
const Epetra_Map &
Ginla::EpetraFVM::KeoPreconditioner::
OperatorDomainMap () const
{
    TEUCHOS_ASSERT( !keoPrec_.is_null() );
    return keoPrec_->OperatorDomainMap();
}
// =============================================================================
const Epetra_Map &
Ginla::EpetraFVM::KeoPreconditioner::
OperatorRangeMap () const
{
    TEUCHOS_ASSERT( !keoPrec_.is_null() );
    return keoPrec_->OperatorRangeMap();
}
// =============================================================================
void
Ginla::EpetraFVM::KeoPreconditioner::
rebuild()
{
    Teuchos::TimeMonitor tm(*rebuildTime_);
    // -------------------------------------------------------------------------
    // rebuild the keo
    keoFactory_->buildKeo( *keoPrec_ );
    keoPrec_->Scale( -1.0 );
    // -------------------------------------------------------------------------
    // regularization
    double mu = keoFactory_->getMvpParameters()->getValue( "mu" );
    if ( fabs( mu ) < 1.0e-5 )
    {
        Teuchos::TimeMonitor tm(*regularizationTime_);
        // Add a regularization to the diagonal.
        Epetra_Vector e( keoPrec_->DomainMap() );
        TEUCHOS_ASSERT_EQUALITY( 0, e.PutScalar( 1.0 ) );

        Epetra_Vector diag( keoPrec_->DomainMap() );
        TEUCHOS_ASSERT_EQUALITY( 0, keoPrec_->ExtractDiagonalCopy( diag ) );
        TEUCHOS_ASSERT_EQUALITY( 0, diag.Update( 1.0e-3, e, 1.0 ) );
        TEUCHOS_ASSERT_EQUALITY( 0, keoPrec_->ReplaceDiagonalValues( diag ) );
    }
    // -------------------------------------------------------------------------
    // Rebuild preconditioner for this object. Not to be mistaken for the
    // object itself.
    switch ( invType_ )
    {
      case Ilu:
          return this->rebuildIlu_();
          break;
      case Ml:
          return this->rebuildMl_();
          break;
      default:
          TEUCHOS_TEST_FOR_EXCEPTION( true,
                                      std::invalid_argument,
                                      "Illegal value of the invType \"" << invType_ << "\"."
                                    );
          break;
    }
    // -------------------------------------------------------------------------

    return;
}
// =============================================================================
void
Ginla::EpetraFVM::KeoPreconditioner::
rebuild( const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams,
         const Teuchos::Tuple<double,3>                  & scaling
       )
{

    keoFactory_->updateParameters( mvpParams, scaling );
    this->rebuild();
    return;
}
// =============================================================================
void
Ginla::EpetraFVM::KeoPreconditioner::
rebuildMl_()
{
    Teuchos::TimeMonitor tm(*rebuildMlTime_);
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
    Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec =
        Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(*keoPrec_, MLList) );
    TEUCHOS_ASSERT( !MLPrec.is_null() );
//     MLPrec->PrintUnused(0);

    // Create the Belos preconditioned operator from the preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    belosPrec_ = Teuchos::rcp( new Belos::EpetraPrecOp( MLPrec ) );

    return;
}
// =============================================================================
void
Ginla::EpetraFVM::KeoPreconditioner::
rebuildIlu_()
{
    Teuchos::TimeMonitor tm(*rebuildIluTime_);
    // set the matrix the linear problem
    keoProblem_->SetOperator( &*keoPrec_ );

    // do the factorizations
    Amesos Factory;
    keoSolver_ = Teuchos::rcp( Factory.Create( "Klu", *keoProblem_ ) );

    // do symbolic and numerical factorizations
    // TODO reuse symbolic factorization
    TEUCHOS_ASSERT_EQUALITY( 0, keoSolver_->SymbolicFactorization() );
    TEUCHOS_ASSERT_EQUALITY( 0, keoSolver_->NumericFactorization() );

    return;
}
// =============================================================================
