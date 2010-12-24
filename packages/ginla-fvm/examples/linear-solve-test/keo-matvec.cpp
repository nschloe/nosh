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

#include <time.h>

// =============================================================================
typedef double                           ST;
typedef Epetra_MultiVector               MV;
typedef Epetra_Operator                  OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
// =============================================================================
timespec
diff(timespec start, timespec end)
{
    timespec temp;
    if ( (end.tv_nsec-start.tv_nsec) < 0 ) {
            temp.tv_sec  = end.tv_sec - start.tv_sec - 1;
            temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    } else {
            temp.tv_sec  = end.tv_sec  - start.tv_sec;
            temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
    return temp;
}
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
      Teuchos::RCP<Ginla::EpetraFVM::StkMesh> mesh = Teuchos::null;

      if ( eComm->MyPID() == 0 )
          std::cout << "Reading..." << std::endl;
      timespec time1, time2;
      int temp;
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);

      Ginla::EpetraFVM::StkMeshRead( *eComm,
                                      inputFileName,
                                      z,
                                      mesh,
                                      problemParameters
                                    );

      clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &time2 );
      if ( eComm->MyPID() == 0 )
          std::cout << "Done. " << diff(time1,time2).tv_sec << "." << setw(9) << setfill('0') << diff(time1,time2).tv_nsec << "s" << std::endl;

      double mu = problemParameters.get<double> ( "mu" );
      mu = 1.0e-3;
      Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp =
              Teuchos::rcp ( new Ginla::MagneticVectorPotential::Z ( mesh, mu ) );

      Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
          Teuchos::rcp( new LOCA::ParameterVector() );
      mvpParameters->addParameter( "mu", mu );

      // create the kinetic energy operator
      Teuchos::RCP<Ginla::EpetraFVM::KeoFactory> keoFactory =
              Teuchos::rcp( new Ginla::EpetraFVM::KeoFactory( mesh, mvp ) );
      Teuchos::RCP<Epetra_FECrsMatrix> keoMatrix =
              Teuchos::rcp( new Epetra_FECrsMatrix( Copy, keoFactory->buildKeoGraph() ) );
      Teuchos::Tuple<double,3> scaling( Teuchos::tuple(1.0,1.0,1.0) );
      keoFactory->buildKeo( *keoMatrix, mvpParameters, scaling );

      // Make sure the matrix is indeed positive definite, and not
      // negative definite. Belos needs that (2010-11-05).
      keoMatrix->Scale( -1.0 );

      // create initial guess and right-hand side
      Teuchos::RCP<Epetra_Vector> epetra_x =
              Teuchos::rcp( new Epetra_Vector( keoMatrix->OperatorDomainMap() ) );
      Teuchos::RCP<Epetra_MultiVector> epetra_b =
              Teuchos::rcp( new Epetra_Vector( keoMatrix->OperatorRangeMap(), 1 ) );
      epetra_b->Random();
      // -----------------------------------------------------------------------
      // perform matrix-vector products

      if ( eComm->MyPID() == 0 )
          std::cout << "MV product..." << std::endl;
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);

      for ( unsigned int k=0; k<2e3; k++ )
          keoMatrix->Apply( *epetra_x, *epetra_b );

      clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &time2 );
      if ( eComm->MyPID() == 0 )
          std::cout << "Done. " << diff(time1,time2).tv_sec << "." << setw(9) << setfill('0') << diff(time1,time2).tv_nsec << "s" << std::endl;

//       clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tS);
//       std::cout << "Done. Time taken is: " << tS.tv_sec << "s, " << tS.tv_nsec << "ns" << std::endl;
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

