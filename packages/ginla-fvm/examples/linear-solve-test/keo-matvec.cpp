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

#include "Ginla_EpetraFVM_StkMeshReader.hpp"

#include "Ginla_EpetraFVM_State.hpp"
#include "Ginla_EpetraFVM_ModelEvaluator.hpp"
#include "Ginla_EpetraFVM_KeoFactory.hpp"
#include "Ginla_IO_StateWriter.hpp"
#include "Ginla_IO_StatsWriter.hpp"
#include "Ginla_IO_NoxObserver.hpp"
#include "Ginla_IO_SaveEigenData.hpp"
#include "Ginla_MagneticVectorPotential_Custom.hpp"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_TimeMonitor.hpp>

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
      if ( eComm->MyPID() == 0 )
          std::cout << "# " << eComm->NumProc() << " processes" << std::endl;
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

      int frequency = 5;
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
      Teuchos::RCP<Epetra_MultiVector>    mvpValues = Teuchos::null;
      Teuchos::RCP<Epetra_Vector>         thickness = Teuchos::null;
      Teuchos::RCP<Ginla::EpetraFVM::StkMesh> mesh = Teuchos::null;

      // read the file
      Teuchos::RCP<Teuchos::Time> readTime = Teuchos::TimeMonitor::getNewTimer("Data I/O");
      {
      Teuchos::TimeMonitor tm(*readTime);
      Ginla::EpetraFVM::StkMeshRead( *eComm,
                                      inputFileName,
                                      z,
                                      mvpValues,
                                      thickness,
                                      mesh,
                                      problemParameters
                                    );
      }

      Teuchos::RCP<Teuchos::Time> mvpConstructTime = Teuchos::TimeMonitor::getNewTimer("MVP construction");
      Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp;
      double mu;
      {
          Teuchos::TimeMonitor tm(*mvpConstructTime);
          mu = problemParameters.get<double> ( "mu" );
          mu = 1.0e-3;
          mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential::Custom ( mesh, mvpValues, mu ) );
          //mvp->initializeEdgeMidpointProjectionCache_();
      }

      Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
          Teuchos::rcp( new LOCA::ParameterVector() );
      mvpParameters->addParameter( "mu", mu );

      Teuchos::RCP<Ginla::EpetraFVM::KeoFactory> keoFactory =
          Teuchos::rcp( new Ginla::EpetraFVM::KeoFactory( mesh, thickness, mvp ) );

      // Precompute FVM entities. Not actually necessary as it's triggered automatically
      // when needed, but for timing purposes put it here.
      Teuchos::RCP<Teuchos::Time> fvmEntitiesConstructTime = Teuchos::TimeMonitor::getNewTimer("FVM entities construction");
      {
          Teuchos::TimeMonitor tm(*fvmEntitiesConstructTime);
          mesh->computeFvmEntities_();
      }

      Teuchos::RCP<Teuchos::Time> graphConstructTime = Teuchos::TimeMonitor::getNewTimer("Graph construction");
      Teuchos::RCP<Epetra_FECrsGraph> keoGraph;
      {
          Teuchos::TimeMonitor tm(*graphConstructTime);
          keoGraph = Teuchos::rcp( new Epetra_FECrsGraph( keoFactory->buildKeoGraph() ) );
      }

      // create the kinetic energy operator
      Teuchos::RCP<Epetra_FECrsMatrix> keoMatrix;
      keoMatrix = Teuchos::rcp( new Epetra_FECrsMatrix( Copy, *keoGraph ) );
      Teuchos::Tuple<double,3> scaling( Teuchos::tuple(1.0,1.0,1.0) );
      Teuchos::RCP<Teuchos::Time> keoConstructTime = Teuchos::TimeMonitor::getNewTimer("Matrix construction");
      {
          Teuchos::TimeMonitor tm(*keoConstructTime);
          keoFactory->updateParameters( mvpParameters, scaling );
          keoFactory->buildKeo( *keoMatrix );
      }
      // Make sure the matrix is indeed positive definite, and not
      // negative definite. Belos needs that (2010-11-05).
      TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->Scale( -1.0 ) );

      // create initial guess and right-hand side
      Teuchos::RCP<Epetra_Vector> epetra_x =
              Teuchos::rcp( new Epetra_Vector( keoMatrix->OperatorDomainMap() ) );
      Teuchos::RCP<Epetra_MultiVector> epetra_b =
              Teuchos::rcp( new Epetra_Vector( keoMatrix->OperatorRangeMap(), 1 ) );
      TEUCHOS_ASSERT_EQUALITY( 0, epetra_x->Random() );

      // -----------------------------------------------------------------------
      // perform several tests
      std::cout << "Process " << eComm->MyPID() << " has " << keoMatrix->NumMyNonzeros() << " nonzeros." << std::endl;

      Teuchos::RCP<Teuchos::Time> mvTime = Teuchos::TimeMonitor::getNewTimer("Matrix-vector multiplication");
      {
          Teuchos::TimeMonitor tm(*mvTime);
          TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->Apply( *epetra_x, *epetra_b ) );
      }

      // print timing data
      Teuchos::TimeMonitor::summarize();
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

    return status==0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================

