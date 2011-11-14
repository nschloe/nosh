// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include <ml_epetra_preconditioner.h>

#include "Ginla_StkMeshReader.hpp"

#include "Ginla_State.hpp"
#include "Ginla_ModelEvaluator.hpp"
#include "Ginla_KeoFactory.hpp"
#include "Ginla_MagneticVectorPotential.hpp"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <EpetraExt_readEpetraLinearSystem.h>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <Amesos_BaseSolver.h>
#include <Amesos.h>

// =============================================================================
int main ( int argc, char *argv[] )
{
  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init( &argc, &argv );
#endif

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm =
      Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm>  eComm =
           Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

    const Teuchos::RCP<Teuchos::FancyOStream> out =
        Teuchos::VerboseObjectBase::getDefaultOStream();

    bool success = true;
    try
    {
      // ===========================================================================
      // handle command line arguments
      Teuchos::CommandLineProcessor My_CLP;

      My_CLP.setDocString (
          "Linear solver testbed for KEO and Jacobian operator.\n"
      );

      std::string inputFileName = "";
      My_CLP.setOption ( "input", &inputFileName,
                         "Input state file", true
                       );

      // print warning for unrecognized arguments
      My_CLP.recogniseAllOptions ( true );

      // do throw exceptions
      My_CLP.throwExceptions ( true );

      // finally, parse the command line
      My_CLP.parse ( argc, argv );
      // =========================================================================
      // Read the data from the file.
      Teuchos::ParameterList data;
      Ginla::StkMeshRead( *eComm, inputFileName, data );

      // Cast the data into something more accessible.
      Teuchos::RCP<Ginla::StkMesh>     & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );
      Teuchos::RCP<Epetra_Vector>      & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );
      Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::ParameterList           & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

      double mu = problemParameters.get<double> ( "mu" );
      Teuchos::RCP<Ginla::MagneticVectorPotential> mvp =
          Teuchos::rcp ( new Ginla::MagneticVectorPotential ( mesh, mvpValues, mu ) );

      Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
          Teuchos::rcp( new LOCA::ParameterVector() );
      mvpParameters->addParameter( "mu", mu );

      // create the kinetic energy operator
      Teuchos::RCP<Ginla::KeoFactory> keoFactory =
              Teuchos::rcp( new Ginla::KeoFactory( mesh, thickness, mvp ) );
      Epetra_FECrsMatrix keoMatrix( Copy, keoFactory->buildKeoGraph() );
      keoFactory->updateParameters( mvpParameters );
      keoFactory->buildKeo( keoMatrix );

      // create initial guess and right-hand side
      Teuchos::RCP<Epetra_Vector> epetra_x = Teuchos::rcp( new Epetra_Vector( keoMatrix.OperatorDomainMap() ) );
      Teuchos::RCP<Epetra_Vector> epetra_b = Teuchos::rcp( new Epetra_Vector( keoMatrix.OperatorRangeMap() ) );
      epetra_b->Random();

      // build the problem
      Epetra_LinearProblem problem( &keoMatrix, &*epetra_x, &*epetra_b );

      // -----------------------------------------------------------------------
//      // Stratimikos.
//      // Thyra glue
//      Teuchos::RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp( keo );
//      Teuchos::RCP<Thyra::VectorBase<double> > x = Thyra::create_Vector( epetra_x, A->domain() );
//      Teuchos::RCP<const Thyra::VectorBase<double> > b = Thyra::create_Vector( epetra_b, A->range() );
//
//      Teuchos::RCP<Teuchos::FancyOStream>
//        out = Teuchos::VerboseObjectBase::getDefaultOStream();
//
//      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
//
//      std::cout << "a" << std::endl;
//
//      // read parameters from file
//      Teuchos::updateParametersFromXmlFile( "./stratimikos.xml", &*linearSolverBuilder.getNonconstParameterList() );
//
//      std::cout << "b" << std::endl;
//
//      // Create a linear solver factory given information read from the
//      // parameter list.
//      Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
//        linearSolverBuilder.createLinearSolveStrategy("");
//
//      // Setup output stream and the verbosity level
//      lowsFactory->setOStream( out );
//      lowsFactory->setVerbLevel( Teuchos::VERB_LOW );
//
//      // Create a linear solver based on the forward operator A
//      Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows =
//        Thyra::linearOpWithSolve(*lowsFactory, A);
//
//      // Solve the linear system (note: the initial guess in 'x' is critical)
//      Thyra::SolveStatus<double> status =
//        Thyra::solve<double>(*lows, Thyra::NOTRANS, *b, x.ptr());
//      *out << "\nSolve status:\n" << status;

      // -----------------------------------------------------------------------
      // build the AztecOO solver
      AztecOO solver( problem );
      // make sure the problem is symmetric
      problem.AssertSymmetric();

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // ML part
      Teuchos::ParameterList MLList;
      ML_Epetra::SetDefaults( "SA", MLList );
      MLList.set("ML output", 10);
      MLList.set("max levels", 10);
      MLList.set("increasing or decreasing", "increasing");
      MLList.set("aggregation: type", "Uncoupled");
      MLList.set("smoother: type", "Chebyshev");
      MLList.set("smoother: sweeps", 3);
      MLList.set("smoother: pre or post", "both");
      MLList.set("coarse: type", "Amesos-KLU");
      MLList.set("PDE equations", 2);
      Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec =
                  Teuchos::rcp( new ML_Epetra::MultiLevelPreconditioner(keoMatrix, MLList) );
      MLPrec->PrintUnused(0);
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      solver.SetAztecOption(AZ_precond, AZ_none);
      //solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
      //solver.SetPrecOperator( MLPrec.getRawPtr() );
      //solver.SetAztecOption(AZ_solver, AZ_gmres);
      solver.SetAztecOption(AZ_solver, AZ_bicgstab);
      //solver.SetAztecOption(AZ_solver, AZ_cg);
      //solver.SetAztecOption(AZ_scaling, 8);
      //solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
      //solver.SetAztecOption(AZ_output, 1);
      //solver.SetAztecOption(AZ_reorder, 0);
      //solver.SetAztecOption(AZ_graph_fill, 3);
      //solver.SetAztecOption(AZ_overlap, 0);
      //solver.SetAztecOption(AZ_poly_ord, 9);
      //solver.SetAztecParam(AZ_ilut_fill, 4.0);
      //solver.SetAztecParam(AZ_drop, 0.0);
      solver.SetAztecOption(AZ_output, 5);
      //double rthresh = 1.4;
      //cout << "Rel threshold = " << rthresh << endl;
      //solver.SetAztecParam(AZ_rthresh, rthresh);
      //double athresh = 10.0;
      //cout << "Abs threshold = " << athresh << endl;
      //solver.SetAztecParam(AZ_athresh, athresh);
      //solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);

      int Niters = 10000;
      //solver.SetAztecOption(AZ_kspace, Niters);

      // do the iteration
      solver.Iterate(Niters, 1.0e-10);

      // compute the residual
      Epetra_Vector bcomp( keoMatrix.OperatorRangeMap() );
      TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix.Apply(*epetra_x, bcomp) );

      Epetra_Vector resid( keoMatrix.OperatorRangeMap() );
      TEUCHOS_ASSERT_EQUALITY( 0, resid.Update(1.0, *epetra_b, -1.0, bcomp, 0.0 ) );

      double residual;
      TEUCHOS_ASSERT_EQUALITY( 0, resid.Norm2(&residual) );
      *out << "Residual    = " << residual << "\n\n" << endl;
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

      *out << "Starting symbolic factorization..." << std::endl;
      Solver->SymbolicFactorization();
      *out << "Starting numeric factorization..." << std::endl;
      Solver->NumericFactorization();
      *out << "Starting solution phase..." << std::endl;
      // solve!
      int ierr = Solver->Solve();
      success = ierr==0;
      // -----------------------------------------------------------------------
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

#ifdef HAVE_MPI
      MPI_Finalize();
#endif

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
