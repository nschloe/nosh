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

#include "Ginla_StkMeshReader.hpp"

#include "Ginla_State.hpp"
#include "Ginla_ModelEvaluator.hpp"
#include "Ginla_KeoContainer.hpp"
#include "Ginla_MagneticVectorPotential_ExplicitValues.hpp"

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

    const Teuchos::RCP<Teuchos::FancyOStream> out =
        Teuchos::VerboseObjectBase::getDefaultOStream();

    bool success = true;
    try
    {
      // ===========================================================================
      *out << "# " << eComm->NumProc() << " processes" << std::endl;
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
      Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );
      Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::ParameterList           & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

      Teuchos::RCP<Teuchos::Time> mvpConstructTime = Teuchos::TimeMonitor::getNewTimer("MVP construction");
      Teuchos::RCP<Ginla::MagneticVectorPotential::ExplicitValues> mvp;
      double mu;
      {
          Teuchos::TimeMonitor tm(*mvpConstructTime);
          mu = problemParameters.get<double> ( "mu" );
          mu = 1.0e-3;
          mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential::ExplicitValues ( mesh, mvpValues, mu ) );
          //mvp->initializeEdgeMidpointProjectionCache_();
      }

      Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
          Teuchos::rcp( new LOCA::ParameterVector() );
      mvpParameters->addParameter( "mu", mu );

      Teuchos::RCP<Ginla::KeoContainer> keoContainer =
          Teuchos::rcp( new Ginla::KeoContainer( mesh, thickness, mvp ) );

      // create the kinetic energy operator
      Teuchos::RCP<Epetra_CrsMatrix> keoMatrix;
      Teuchos::RCP<Teuchos::Time> keoConstructTime = Teuchos::TimeMonitor::getNewTimer("Matrix construction");
      {
          Teuchos::TimeMonitor tm(*keoConstructTime);
          keoContainer->updateParameters( mvpParameters );
          keoMatrix = Teuchos::rcp( new Epetra_CrsMatrix( *keoContainer->getKeo() ) );
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
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

#ifdef HAVE_MPI
      MPI_Finalize();
#endif

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================

