// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Epetra_LinearProblem.h>

#include <ml_epetra_preconditioner.h>

//#include "BelosConfigDefs.hpp"
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosMinresSolMgr.hpp>

#include "Ginla_StkMeshReader.hpp"

#include "Ginla_State.hpp"
#include "Ginla_ModelEvaluator.hpp"
#include "Ginla_KeoContainer.hpp"
#include "Ginla_KeoRegularized.hpp"
#include "Ginla_StatsWriter.hpp"
#include "Ginla_NoxObserver.hpp"
#include "Ginla_SaveEigenData.hpp"
#include "Ginla_MagneticVectorPotential_ExplicitValues.hpp"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// =============================================================================
typedef double                           ST;
typedef Epetra_MultiVector               MV;
typedef Epetra_Operator                  OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

enum Operator { JAC, KEO, KEOREG, POISSON1D };
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

      std::string inputFileName( "" );
      My_CLP.setOption("input", &inputFileName, "Input state file", true);

      std::string action( "solve" );
      My_CLP.setOption("action", &action, "Which action to perform with the operator (solve, matvec)");

      std::string solver( "cg" );
      My_CLP.setOption("solver", &solver, "Krylov subspace method (cg, minres, gmres)");

      std::string op( "jac" );
      My_CLP.setOption("operator", &op, "Operator (jac, keo, keoreg, poisson1d)");
//       Operator op = JAC;
//       Operator allOpts[] = {JAC, KEO, KEOREG, POISSON1D};
//       std::string allOptNames[] = {"jac", "keo", "keoreg", "poisson1d"};
//       My_CLP.setOption("operator", &op, 4, allOpts, allOptNames);

      bool verbose = true;
      My_CLP.setOption("verbose", "quiet", &verbose, "Print messages and results.");

      bool usePrec = true;
      My_CLP.setOption("prec", "noprec", &usePrec, "Use a preconditioner.");

      int frequency = 10;
      My_CLP.setOption("frequency", &frequency, "Solvers frequency for printing residuals (#iters).");

      // print warning for unrecognized arguments
      My_CLP.recogniseAllOptions ( true );

      // finally, parse the command line
      TEUCHOS_ASSERT_EQUALITY( My_CLP.parse ( argc, argv ),
                               Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL
                             );
      // =========================================================================
      // Read the data from the file.
      Teuchos::ParameterList data;
      Teuchos::RCP<Teuchos::Time> readTime = Teuchos::TimeMonitor::getNewTimer("Data I/O");
      {
      Teuchos::TimeMonitor tm(*readTime);
      Ginla::StkMeshRead( *eComm, inputFileName, data );
      }

      // Cast the data into something more accessible.
      Teuchos::RCP<Ginla::StkMesh>     & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );
      Teuchos::RCP<Epetra_Vector>      & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::RCP<const Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<const Epetra_MultiVector>() );
      Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::ParameterList           & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

      Teuchos::RCP<Ginla::MagneticVectorPotential::ExplicitValues> mvp;
      double mu;
      Teuchos::RCP<Teuchos::Time> mvpConstructTime = Teuchos::TimeMonitor::getNewTimer("MVP construction");
      {
          Teuchos::TimeMonitor tm(*mvpConstructTime);
//          mu = problemParameters.get<double> ( "mu" );
          mu = 1.0;
          mvp = Teuchos::rcp(new Ginla::MagneticVectorPotential::ExplicitValues(mesh, mvpValues, mu));
      }

      Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
          Teuchos::rcp( new LOCA::ParameterVector() );
      mvpParameters->addParameter( "mu", mu );

      Teuchos::RCP<Teuchos::Time> keoContainerConstructTime =
          Teuchos::TimeMonitor::getNewTimer("Keo container construction");
      Teuchos::RCP<Ginla::KeoContainer> keoContainer;
      {
          Teuchos::TimeMonitor tm(*keoContainerConstructTime);
          keoContainer = Teuchos::rcp( new Ginla::KeoContainer( mesh, thickness, mvp ) );
      }

      // create Jacobian
      Teuchos::RCP<Teuchos::Time> operatorConstructTime =
          Teuchos::TimeMonitor::getNewTimer("Operator construction");
      Teuchos::RCP<const Epetra_Operator> A;
      {
          Teuchos::TimeMonitor tm(*operatorConstructTime);
          if (op.compare("jac") == 0)
              A = Teuchos::rcp(new Ginla::JacobianOperator(mesh, thickness, keoContainer, z));
          else if (op.compare("keo") == 0)
              A = keoContainer->getKeo();
          else if (op.compare("keoreg") == 0)
              A = Teuchos::rcp(new Ginla::KeoRegularized(mesh, thickness, keoContainer, z));
          else if (op.compare("poisson1d") == 0)
          {
              int n = 1000;
              // Build the matrix (-1,2,-1).
              const double neg_one = -1.0;
              const double two = 2.0;
              Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(n, 0, *eComm));
              Teuchos::RCP<Epetra_CrsMatrix> AMatrix =
                  Teuchos::rcp(new Epetra_CrsMatrix(Copy, *map, 3));
              for (int k=0; k < map->NumMyElements(); k++)
              {
                int k_min1 = k - 1;
                int k_plu1 = k + 1;
                if (k>0)
                  AMatrix->InsertMyValues(k, 1, &neg_one, &k_min1);
                AMatrix->InsertMyValues(k, 1, &two, &k);
                if (k<n-1)
                  AMatrix->InsertMyValues(k, 1, &neg_one, &k_plu1);
              }
              AMatrix->FillComplete();
              A = AMatrix;
          }
          else
              TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Unknown operator \"" << op << "\"." );
      }

      // create initial guess and right-hand side
      Teuchos::RCP<Epetra_Vector> epetra_x =
              Teuchos::rcp( new Epetra_Vector( A->OperatorDomainMap() ) );
      Teuchos::RCP<Epetra_MultiVector> epetra_b =
              Teuchos::rcp( new Epetra_Vector( A->OperatorRangeMap() ) );
      // epetra_b->Random();
      epetra_b->PutScalar( 1.0 );


      if (action.compare("matvec") == 0)
      {
        epetra_x->PutScalar( 1.0 );
        Teuchos::RCP<Teuchos::Time> mvTime = Teuchos::TimeMonitor::getNewTimer("Operator apply");
        {
          Teuchos::TimeMonitor tm(*mvTime);
          TEUCHOS_ASSERT_EQUALITY(0, A->Apply(*epetra_x, *epetra_b));
        }
        // print timing data
        Teuchos::TimeMonitor::summarize();
      }
      else if (action.compare("solve") == 0)
      {
        // -----------------------------------------------------------------------
        // Belos part
        Teuchos::ParameterList belosList;
        // Relative convergence tolerance requested
        belosList.set( "Convergence Tolerance", 1.0e-12 );
        if (verbose) {
          belosList.set( "Verbosity",
                        Belos::Errors +
                        Belos::Warnings +
                        Belos::IterationDetails +
                        Belos::FinalSummary +
                        Belos::Debug +
                        Belos::TimingDetails +
                        Belos::StatusTestDetails
                      );
          if (frequency > 0)
            belosList.set("Output Frequency", frequency);
        }
        else
          belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

        belosList.set( "Output Style", (int)Belos::Brief ); // Belos::General, Belos::Brief
        belosList.set( "Maximum Iterations", 1000 );

        // Construct an unpreconditioned linear problem instance.
        Belos::LinearProblem<double,MV,OP> problem( A, epetra_x, epetra_b );
        bool set = problem.setProblem();
        TEUCHOS_TEST_FOR_EXCEPTION(!set,
                                  std::logic_error,
                                  "ERROR:  Belos::LinearProblem failed to set up correctly!" );
        // -----------------------------------------------------------------------
        // create preconditioner
        Teuchos::RCP<Teuchos::Time> precConstructTime = Teuchos::TimeMonitor::getNewTimer("Prec construction");
        if ( usePrec )
        {
            Teuchos::TimeMonitor tm(*precConstructTime);

            // create the preconditioner
            Teuchos::RCP<Ginla::KeoRegularized> keoReg =
                Teuchos::rcp(new Ginla::KeoRegularized(mesh, thickness, keoContainer, z));

            // initialize its (approximate) inverse
            keoReg->rebuildInverse();

            // Create the Belos preconditioned operator from the preconditioner.
            // NOTE:  This is necessary because Belos expects an operator to apply the
            //        preconditioner with Apply() NOT ApplyInverse().
            Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp(new Belos::EpetraPrecOp(keoReg));
            problem.setLeftPrec( belosPrec );
        }
        // -----------------------------------------------------------------------
        // Create an iterative solver manager.
        Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver;
        if (solver.compare("cg") == 0)
        {
          belosList.set( "Assert Positive Definiteness", false );
          newSolver =
            Teuchos::rcp(new Belos::PseudoBlockCGSolMgr<double,MV,OP>(Teuchos::rcp(&problem,false),
                                                                      Teuchos::rcp(&belosList,false)
                                                                      )
                        );
        }
        else if (solver.compare("minres") == 0)
        {
          newSolver =
            Teuchos::rcp(new Belos::MinresSolMgr<double,MV,OP>(Teuchos::rcp(&problem,false),
                                                              Teuchos::rcp(&belosList,false)
                                                              )
                        );
        }
        else if (solver.compare("gmres") == 0)
        {
          newSolver =
            Teuchos::rcp(new Belos::PseudoBlockGmresSolMgr<double,MV,OP>(Teuchos::rcp(&problem,false),
                                                                        Teuchos::rcp(&belosList,false)
                                                                        )
                        );
        }
        else
        {
          TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Unknown solver type \"" << solver << "\"." );
        }

        // Perform solve
        Teuchos::RCP<Teuchos::Time> solveTime =
          Teuchos::TimeMonitor::getNewTimer("Jacobian solve");
        {
            Teuchos::TimeMonitor tm(*solveTime);
            Belos::ReturnType ret = newSolver->solve();
            success = ret==Belos::Converged;
        }

//         *out << newSolver->getNumIters() << std::endl;
//         // Compute actual residuals.
//         bool badRes = false;
//         Teuchos::Array<double> actual_resids( 1 );
//         Teuchos::Array<double> rhs_norm( 1 );
//         Epetra_Vector resid( keoMatrix->OperatorRangeMap() );
//         OPT::Apply( *keoMatrix, *epetra_x, resid );
//         MVT::MvAddMv( -1.0, resid, 1.0, *epetra_b, resid );
//         MVT::MvNorm( resid, actual_resids );
//         MVT::MvNorm( *epetra_b, rhs_norm );
//         if (proc_verbose) {
//           std::cout<< "---------- Actual Residuals (normalized) ----------" <<std::endl<<std::endl;
//           for ( int i=0; i<1; i++) {
//             double actRes = actual_resids[i]/rhs_norm[i];
//             std::cout << "Problem " << i << " : \t" << actRes << std::endl;
//             if (actRes > 1.0e-10) badRes = true;
//           }
//         }
      }
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

#ifdef HAVE_MPI
      MPI_Finalize();
#endif

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
