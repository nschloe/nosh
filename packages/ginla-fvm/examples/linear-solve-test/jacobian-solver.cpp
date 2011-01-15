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

//#include "BelosConfigDefs.hpp"
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>

#include <boost/filesystem.hpp>

#include "Ginla_EpetraFVM_StkMeshReader.h"

#include "Ginla_EpetraFVM_State.h"
#include "Ginla_EpetraFVM_ModelEvaluator.h"
#include "Ginla_EpetraFVM_KeoFactory.h"
#include "Ginla_EpetraFVM_KeoPreconditioner.h"
#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_StatsWriter.h"
#include "Ginla_IO_NoxObserver.h"
#include "Ginla_IO_SaveEigenData.h"
#include "Ginla_MagneticVectorPotential_Spherical.h"
#include "Ginla_MagneticVectorPotential_Z.h"
#include "Ginla_MagneticVectorPotential_MagneticDot.h"

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

    int status = 0;
    Belos::ReturnType ret;

    try
    {
      // ===========================================================================
      // handle command line arguments
      Teuchos::CommandLineProcessor My_CLP;

      My_CLP.setDocString (
              "Linear solver testbed for KEO and Jacobian operator.\n"
      );

      std::string inputFileName( "" );
      My_CLP.setOption ( "input", &inputFileName, "Input state file", true );

      bool verbose = true;
      My_CLP.setOption("verbose","quiet",&verbose,"Print messages and results.");

      bool isPrec = true;
      My_CLP.setOption("prec","noprec",&isPrec,"Use a preconditioner.");

      int frequency = 10;
      My_CLP.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");

      // print warning for unrecognized arguments
      My_CLP.recogniseAllOptions ( true );

      // finally, parse the command line
      TEUCHOS_ASSERT_EQUALITY( My_CLP.parse ( argc, argv ),
                               Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL
                             );
      // =========================================================================
      Teuchos::ParameterList              problemParameters;
      Teuchos::RCP<Epetra_Vector>         z = Teuchos::null;
      Teuchos::RCP<Ginla::EpetraFVM::StkMesh> mesh = Teuchos::null;

      if ( eComm->MyPID() == 0 )
          std::cout << "Reading..." << std::endl;

      Teuchos::RCP<Teuchos::Time> readTime = Teuchos::TimeMonitor::getNewTimer("Data I/O");
      {
      Teuchos::TimeMonitor tm(*readTime);
      Ginla::EpetraFVM::StkMeshRead( *eComm,
                                      inputFileName,
                                      z,
                                      mesh,
                                      problemParameters
                                    );
      }

      Teuchos::RCP<Ginla::MagneticVectorPotential::Z> mvp;
      double mu;
      Teuchos::RCP<Teuchos::Time> mvpConstructTime = Teuchos::TimeMonitor::getNewTimer("MVP construction");
      {
          Teuchos::TimeMonitor tm(*mvpConstructTime);
          mu = problemParameters.get<double> ( "mu" );
          mu = 1.0e-2;
          mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential::Z ( mesh, mu ) );
          mvp->initializeEdgeMidpointProjectionCache_();
      }

      Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
          Teuchos::rcp( new LOCA::ParameterVector() );
      mvpParameters->addParameter( "mu", mu );

      // Precompute FVM entities. Not actually necessary as it's triggered automatically
      // when needed, but for timing purposes put it here.
      Teuchos::RCP<Teuchos::Time> fvmEntitiesConstructTime = Teuchos::TimeMonitor::getNewTimer("FVM entities construction");
      {
          Teuchos::TimeMonitor tm(*fvmEntitiesConstructTime);
          mesh->computeFvmEntities_();
      }

      Teuchos::RCP<Ginla::EpetraFVM::KeoFactory> keoFactory =
          Teuchos::rcp( new Ginla::EpetraFVM::KeoFactory( mesh, mvp ) );

      Teuchos::RCP<Epetra_FECrsGraph> keoGraph;
      Teuchos::RCP<Teuchos::Time> graphConstructTime = Teuchos::TimeMonitor::getNewTimer("Graph construction");
      {
          Teuchos::TimeMonitor tm(*graphConstructTime);
          keoGraph = Teuchos::rcp( new Epetra_FECrsGraph( keoFactory->buildKeoGraph() ) );
      }

      // create Jacobian
      Teuchos::RCP<Teuchos::Time> jacobianConstructTime = Teuchos::TimeMonitor::getNewTimer("Jacobian construction");
      Teuchos::RCP<Ginla::EpetraFVM::JacobianOperator> jac;
      {
          Teuchos::TimeMonitor tm(*jacobianConstructTime);
          // create the jacobian operator
          jac = Teuchos::rcp( new Ginla::EpetraFVM::JacobianOperator( mesh, mvp, z ) );
      }

      // create initial guess and right-hand side
      Teuchos::RCP<Epetra_Vector> epetra_x =
              Teuchos::rcp( new Epetra_Vector( jac->OperatorDomainMap() ) );
      Teuchos::RCP<Epetra_MultiVector> epetra_b =
              Teuchos::rcp( new Epetra_Vector( jac->OperatorRangeMap(), 1 ) );
      epetra_b->Random();


      // -----------------------------------------------------------------------
      // Belos part
      ParameterList belosList;
      belosList.set( "Convergence Tolerance", 1.0e-15 );  // Relative convergence tolerance requested
      if (verbose) {
        belosList.set( "Verbosity",
                       Belos::Errors +
                       Belos::Warnings +
                       Belos::IterationDetails +
                       Belos::FinalSummary +
                       Belos::Debug +
                       Belos::TimingDetails //+
//                       Belos::StatusTestDetails
                     );
        if (frequency > 0)
          belosList.set( "Output Frequency", frequency );
      }
      else
        belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

      // Only print on the zero processor
      const bool proc_verbose = verbose && (eComm->MyPID()==0);

      // Construct an unpreconditioned linear problem instance.
      Belos::LinearProblem<double,MV,OP> problem( jac, epetra_x, epetra_b );
      bool set = problem.setProblem();
      if (set == false)
      {
          if (proc_verbose)
              std::cout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
          return -1;
      }
      // -----------------------------------------------------------------------
      // create preconditioner
      // create preconditioner
      Teuchos::RCP<Teuchos::Time> precConstructTime = Teuchos::TimeMonitor::getNewTimer("Prec construction");
      if ( isPrec )
      {
          Teuchos::TimeMonitor tm(*precConstructTime);
          // create the jacobian operator
          Teuchos::RCP<Ginla::EpetraFVM::KeoPreconditioner> prec =
              Teuchos::rcp( new Ginla::EpetraFVM::KeoPreconditioner( mesh, mvp ) );

          // Create the Belos preconditioned operator from the preconditioner.
          // NOTE:  This is necessary because Belos expects an operator to apply the
          //        preconditioner with Apply() NOT ApplyInverse().
          Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( prec ) );
          problem.setLeftPrec( belosPrec );
      }
      // -----------------------------------------------------------------------
      // Create an iterative solver manager.
      Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver
              = Teuchos::rcp( new Belos::PseudoBlockCGSolMgr<double,MV,OP>( Teuchos::rcp(&problem,false),
                                                                            Teuchos::rcp(&belosList,false)
                                                                          )
                   );
//       RCP< Belos::SolverManager<double,MV,OP> > newSolver
//               = rcp( new Belos::PseudoBlockGmresSolMgr<double,MV,OP>( rcp(&problem,false),
//                                                                       rcp(&belosList,false)
//                                                                     )
//                    );

      // Perform solve
      ret = newSolver->solve();
      // -----------------------------------------------------------------------
      //
      // Compute actual residuals.
      //
//       bool badRes = false;
//       std::vector<double> actual_resids( 1 );
//       std::vector<double> rhs_norm( 1 );
//       Epetra_Vector resid( keoMatrix->OperatorRangeMap() );
//       OPT::Apply( *keoMatrix, *epetra_x, resid );
//       MVT::MvAddMv( -1.0, resid, 1.0, *epetra_b, resid );
//       MVT::MvNorm( resid, actual_resids );
//       MVT::MvNorm( *epetra_b, rhs_norm );
//       if (proc_verbose) {
//         std::cout<< "---------- Actual Residuals (normalized) ----------" <<std::endl<<std::endl;
//         for ( int i=0; i<1; i++) {
//           double actRes = actual_resids[i]/rhs_norm[i];
//           std::cout << "Problem " << i << " : \t" << actRes << std::endl;
//           if (actRes > 1.0e-10) badRes = true;
//         }
//       }
      // -----------------------------------------------------------------------
    }
    catch ( std::exception & e )
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << e.what() << std::endl;
        status += 10;
    }
    catch ( std::string & e )
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << e << std::endl;
        status += 10;
    }
    catch ( const char * e )
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << e << std::endl;
        status += 10;
    }
    catch ( int e )
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << "Caught unknown exception code " << e <<  "." << std::endl;
        status += 10;
    }
    catch (...)
    {
        if ( eComm->MyPID() == 0 )
            std::cerr << "Caught unknown exception." << std::endl;
        status += 10;
    }

#ifdef HAVE_MPI
      MPI_Finalize();
#endif

    return ret==Belos::Converged ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
