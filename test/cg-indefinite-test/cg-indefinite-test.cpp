#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <Amesos.h>

#include <ml_epetra_preconditioner.h>

#include <boost/filesystem.hpp>

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// =============================================================================
int main ( int argc, char *argv[] )
{
  // Initialize MPI
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_MpiComm> eComm =
    Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
  Teuchos::RCP<Epetra_SerialComm>  eComm =
    Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  bool success = true;
  try
  {
    // define the map
    int N = 1e1;
    Teuchos::RCP<Epetra_Map> map =
            Teuchos::rcp( new Epetra_Map(N,0,*eComm) );

    // create a diagonal matrix
    Teuchos::RCP<Epetra_CrsMatrix> A =
            Teuchos::rcp( new Epetra_CrsMatrix( Copy, *map, 1 ));
    Epetra_Vector diag(*map);
    diag.Random();
    for ( int k=0; k<map->NumMyElements(); k++ )
    {
        double val = 0.5 * diag[k] + 1.0;
        TEUCHOS_ASSERT_EQUALITY( 0, A->InsertGlobalValues( k, 1, &val, &k ) );
    }
    TEUCHOS_ASSERT_EQUALITY( 0, A->FillComplete() );
    int k = 0;
    double val = -1.0;
    TEUCHOS_ASSERT_EQUALITY( 0, A->ReplaceGlobalValues( k, 1, &val, &k ) );

    A->Print( std::cout );

    // create initial guess and right-hand side
    Teuchos::RCP<Epetra_Vector> epetra_x = Teuchos::rcp( new Epetra_Vector( A->DomainMap() ) );
    Teuchos::RCP<Epetra_Vector> epetra_b = Teuchos::rcp( new Epetra_Vector( A->RangeMap() ) );
    epetra_b->Random();
    epetra_b->PutScalar(1.0);

    // build the problem
    Epetra_LinearProblem problem( &*A, &*epetra_x, &*epetra_b );

    // -----------------------------------------------------------------------
    //// Stratimikos.
    //// Thyra glue
    //Teuchos::RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp( keo );
    //Teuchos::RCP<Thyra::VectorBase<double> > x = Thyra::create_Vector( epetra_x, A->domain() );
    //Teuchos::RCP<const Thyra::VectorBase<double> > b = Thyra::create_Vector( epetra_b, A->range() );

    //Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    //std::cout << "a" << std::endl;

    //// read parameters from file
    //Teuchos::updateParametersFromXmlFile( "./stratimikos.xml", &*linearSolverBuilder.getNonconstParameterList() );

    //std::cout << "b" << std::endl;

    //// Create a linear solver factory given information read from the
    //// parameter list.
    //Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
    //  linearSolverBuilder.createLinearSolveStrategy("");

    //// Setup output stream and the verbosity level
    //lowsFactory->setOStream( out );
    //lowsFactory->setVerbLevel( Teuchos::VERB_LOW );

    //// Create a linear solver based on the forward operator A
    //Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows =
    //  Thyra::linearOpWithSolve(*lowsFactory, A);

    //// Solve the linear system (note: the initial guess in 'x' is critical)
    //Thyra::SolveStatus<double> status =
    //  Thyra::solve<double>(*lows, Thyra::NOTRANS, *b, x.ptr());
    //*out << "\nSolve status:\n" << status;

    // -----------------------------------------------------------------------
    // build the AztecOO solver
    AztecOO solver( problem );
    // make sure the problem is symmetric
    problem.AssertSymmetric();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // ML part
//      Teuchos::ParameterList MLList;
//      ML_Epetra::SetDefaults( "SA", MLList );
//      MLList.set("ML output", 10);
//      MLList.set("max levels", 10);
//      MLList.set("increasing or decreasing", "increasing");
//      MLList.set("aggregation: type", "Uncoupled");
//      MLList.set("smoother: type", "Chebyshev");
//      MLList.set("smoother: sweeps", 3);
//      MLList.set("smoother: pre or post", "both");
//      MLList.set("coarse: type", "Amesos-KLU");
//      Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec =
//                  Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(*keoMatrix, MLList) );
//      MLPrec->PrintUnused(0);
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    solver.SetAztecOption(AZ_precond, AZ_none);
    //solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    //solver.SetPrecOperator( MLPrec.getRawPtr() );
    //solver.SetAztecOption(AZ_solver, AZ_gmres);
    solver.SetAztecOption(AZ_solver, AZ_cg);
    //solver.SetAztecOption(AZ_scaling, 8);
    //solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    solver.SetAztecOption(AZ_output, 1);
    //solver.SetAztecOption(AZ_reorder, 0);
    //solver.SetAztecOption(AZ_graph_fill, 3);
    //solver.SetAztecOption(AZ_overlap, 0);
    //solver.SetAztecOption(AZ_poly_ord, 9);
    //solver.SetAztecParam(AZ_ilut_fill, 4.0);
    //solver.SetAztecParam(AZ_drop, 0.0);
    solver.SetAztecOption(AZ_output, 32);
    //double rthresh = 1.4;
    //cout << "Rel threshold = " << rthresh << endl;
    //solver.SetAztecParam(AZ_rthresh, rthresh);
    //double athresh = 10.0;
    //cout << "Abs threshold = " << athresh << endl;
    //solver.SetAztecParam(AZ_athresh, athresh);
    //solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);

    int Niters = N;
    //solver.SetAztecOption(AZ_kspace, Niters);

    // do the iteration
    solver.Iterate(Niters, 1.0e-10);

    // compute the residual
    Epetra_Vector bcomp( A->RangeMap() );
    TEUCHOS_ASSERT_EQUALITY( 0, A->Apply(*epetra_x, bcomp) );

    Epetra_Vector resid( A->RangeMap() );
    TEUCHOS_ASSERT_EQUALITY( 0, resid.Update(1.0, *epetra_b, -1.0, bcomp, 0.0 ) );

    double residual;
    TEUCHOS_ASSERT_EQUALITY( 0, resid.Norm2(&residual) );
    if (eComm->MyPID()==0) cout << "Residual    = " << residual << "\n\n" << endl;
    // -----------------------------------------------------------------------
    // direct solver
    Amesos Factory;
    std::string SolverType = "Klu";
    Teuchos::RCP<Amesos_BaseSolver> Solver =
            Teuchos::rcp( Factory.Create( SolverType, problem ) );
    TEUCHOS_ASSERT( !Solver.is_null() );

    Teuchos::ParameterList List;
    List.set("PrintTiming", true);
    List.set("PrintStatus", true);
    Solver->SetParameters(List);

    if (eComm->MyPID() == 0)
      std::cout << "Starting symbolic factorization..." << std::endl;
    Solver->SymbolicFactorization();
    if (eComm->MyPID() == 0)
      std::cout << "Starting numeric factorization..." << std::endl;
    Solver->NumericFactorization();
    if (eComm->MyPID() == 0)
      std::cout << "Starting solution phase..." << std::endl;
    // solve!
    Solver->Solve();
    // -----------------------------------------------------------------------
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
