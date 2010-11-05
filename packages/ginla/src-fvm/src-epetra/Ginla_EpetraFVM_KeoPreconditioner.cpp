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

#include "Ginla_EpetraFVM_KeoPreconditioner.h"

#include "Ginla_EpetraFVM_KineticEnergyOperator.h"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <ml_epetra_preconditioner.h>

// =============================================================================
// some typdefs for Belos
typedef double                           ST;
typedef Epetra_MultiVector               MV;
typedef Epetra_Operator                  OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
// =============================================================================
Ginla::EpetraFVM::KeoPreconditioner::
KeoPreconditioner( const Teuchos::RCP<VIO::EpetraMesh::Mesh>                   & mesh,
                   const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp
                 ):
        useTranspose_ ( false ),
        comm_( Teuchos::rcpFromRef(mesh->getNodesMap()->Comm()) ),
        keo_( Teuchos::rcp( new Ginla::EpetraFVM::KineticEnergyOperator(mesh, mvp) ) ),
        isKeoIllConditioned_( false )
//        keoProblem_( Teuchos::rcp( new Epetra_LinearProblem() ) ),
//        keoSolver_( Teuchos::null ),
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
    TEUCHOS_ASSERT_EQUALITY( 0, keo_->Apply( X, Y ) );

    if ( isKeoIllConditioned_ )
    {
        // Add a regularization to the diagonal.
        // TODO Put this in the matrix directly.
        Epetra_Vector e( X.Map() );
        e.PutScalar( 1.0e-3 );
        Y.Update( -1.0, e, 1.0 );
    }

    return 0;
}
// =============================================================================
int
Ginla::EpetraFVM::KeoPreconditioner::
ApplyInverse ( const Epetra_MultiVector & X,
                     Epetra_MultiVector & Y
             ) const
{
    bool verbose = true;
    int frequency = 10;
    bool proc_verbose = verbose && (X.Comm().MyPID()==0);

    // Belos part
    ParameterList belosList;
    belosList.set( "Convergence Tolerance", 1.0e-10 );         // Relative convergence tolerance requested
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
                     Belos::TimingDetails + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

    // Construct an unpreconditioned linear problem instance.
    Teuchos::RCP<const Epetra_MultiVector> Xptr = Teuchos::rcpFromRef( X );
    Teuchos::RCP<Epetra_MultiVector> Yptr = Teuchos::rcpFromRef( Y );
    Belos::LinearProblem<double,MV,OP> problem( keo_, Yptr, Xptr );
    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cerr << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    // -------------------------------------------------------------------------
    // create preconditioner
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults( "SA", MLList );
    MLList.set("ML output", 0);
    MLList.set("max levels", 10);
    MLList.set("increasing or decreasing", "increasing");
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("smoother: type", "Chebyshev"); // "block Gauss-Seidel" "Chebyshev"
//      MLList.set("aggregation: threshold", 0.0);
    MLList.set("smoother: sweeps", 3);
    MLList.set("smoother: pre or post", "both");
    MLList.set("coarse: type", "Amesos-KLU");
    MLList.set("PDE equations", 2);
    Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec =
                Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(*keo_->getMatrix(), MLList) );
    MLPrec->PrintUnused(0);

    // create preconditioner object
    Teuchos::RCP<Epetra_Operator> Prec =
            Teuchos::rcp(  new ML_Epetra::MultiLevelPreconditioner(*keo_->getMatrix(), MLList) );
    TEUCHOS_ASSERT( !Prec.is_null() );

    // Create the Belos preconditioned operator from the preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec =
            Teuchos::rcp( new Belos::EpetraPrecOp( Prec ) );
    problem.setLeftPrec( belosPrec );
    // -------------------------------------------------------------------------
    // Create an iterative solver manager.
    Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver =
            Teuchos::rcp( new Belos::BlockCGSolMgr<double,MV,OP>( rcp(&problem,false),
                                                                  rcp(&belosList,false)
                                                                )
                        );

    // Perform solve
    Belos::ReturnType ret = newSolver->solve();

    return ret==Belos::Converged ? 0 : -1;
}
// =============================================================================
//int
//Ginla::EpetraFVM::KeoPreconditioner::
//ApplyInverse2 ( const Epetra_MultiVector & X,
//                     Epetra_MultiVector & Y
//             ) const
//{
//    if ( keoSolver_.is_null() )
//    {
//        // set the matrix the linear problem
//        keoProblem_->SetOperator( &*keo_ );
//
//        // do the factorizations
//        Amesos Factory;
//        keoSolver_ = Teuchos::rcp( Factory.Create( "Klu", *keoProblem_ ) );
//
//        // do symbolic and numerical factorizations
//        // TODO reuse symbolic factorization
//        keoSolver_->SymbolicFactorization();
//        keoSolver_->NumericFactorization();
//    }
//
//    // set left- and right-hand side
//    TEUCHOS_ASSERT( !keoProblem_.is_null() );
//    keoProblem_->SetLHS( &Y );
//    Epetra_MultiVector T = X;
//    keoProblem_->SetRHS( &T );
//
//    // solve and return error code
//    return keoSolver_->Solve();
//}
// =============================================================================
double
Ginla::EpetraFVM::KeoPreconditioner::
NormInf () const
{
    TEST_FOR_EXCEPTION( true,
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
    TEUCHOS_ASSERT( !keo_.is_null() );
    return keo_->OperatorDomainMap();
}
// =============================================================================
const Epetra_Map &
Ginla::EpetraFVM::KeoPreconditioner::
OperatorRangeMap () const
{
    TEUCHOS_ASSERT( !keo_.is_null() );
    return keo_->OperatorRangeMap();
}
// =============================================================================
void
Ginla::EpetraFVM::KeoPreconditioner::
rebuild( const double mu,
         const Teuchos::Tuple<double,3> & scaling
       )
{
    // rebuild the keo
    keo_->setParameters( mu, scaling );
    keo_->assembleMatrix();

    if ( fabs(mu) < 1.0e-5 )
        isKeoIllConditioned_ = true;
    else
        isKeoIllConditioned_ = false;

    return;
}
// =============================================================================
