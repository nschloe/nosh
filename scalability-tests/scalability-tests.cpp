// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Epetra_Vector.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_TimeMonitor.hpp>

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

    try
    {
      // =======================================================================
      if ( eComm->MyPID() == 0 )
          std::cout << "# " << eComm->NumProc() << " processes" << std::endl;
      // =======================================================================
      // Create map.
      // Do strong scaling tests, so keep numGlobalElements independent of
      // the number of processes.
      int numGlobalElements = 1e5;
      int indexBase = 0;
      Teuchos::RCP<Epetra_Map> map =
          Teuchos::rcp( new Epetra_Map ( numGlobalElements, indexBase, *eComm ) );
      // =======================================================================
      // tests on one vector
      Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp( new Epetra_Vector( *map ) );
      u->Random();

      Teuchos::RCP<Teuchos::Time> meanValueTime =
          Teuchos::TimeMonitor::getNewTimer("Vector::MeanValue");
      {
          Teuchos::TimeMonitor tm(*meanValueTime);
          double meanVal;
          TEUCHOS_ASSERT_EQUALITY( 0, u->MeanValue( &meanVal ) );
      }

      Teuchos::RCP<Teuchos::Time> maxValueTime =
          Teuchos::TimeMonitor::getNewTimer("Vector::MaxValue");
      {
          Teuchos::TimeMonitor tm(*maxValueTime);
          double maxValue;
          TEUCHOS_ASSERT_EQUALITY( 0, u->MaxValue( &maxValue ) );
      }

      Teuchos::RCP<Teuchos::Time> minValueTime =
          Teuchos::TimeMonitor::getNewTimer("Vector::MinValue");
      {
          Teuchos::TimeMonitor tm(*minValueTime);
          double minValue;
          TEUCHOS_ASSERT_EQUALITY( 0, u->MinValue( &minValue ) );
      }

      Teuchos::RCP<Teuchos::Time> norm1Time =
          Teuchos::TimeMonitor::getNewTimer("Vector::Norm1");
      {
          Teuchos::TimeMonitor tm(*norm1Time);
          double norm1;
          TEUCHOS_ASSERT_EQUALITY( 0, u->Norm1( &norm1 ) );
      }

      Teuchos::RCP<Teuchos::Time> norm2Time =
          Teuchos::TimeMonitor::getNewTimer("Vector::Norm2");
      {
          Teuchos::TimeMonitor tm(*norm2Time);
          double norm2;
          TEUCHOS_ASSERT_EQUALITY( 0, u->Norm2( &norm2 ) );
      }

      Teuchos::RCP<Teuchos::Time> normInfTime =
          Teuchos::TimeMonitor::getNewTimer("Vector::NormInf");
      {
          Teuchos::TimeMonitor tm(*normInfTime);
          double normInf;
          TEUCHOS_ASSERT_EQUALITY( 0, u->NormInf( &normInf ) );
      }

      Teuchos::RCP<Teuchos::Time> scaleTime =
          Teuchos::TimeMonitor::getNewTimer("Vector::Scale");
      {
          Teuchos::TimeMonitor tm(*scaleTime);
          double alpha = 0.5;
          TEUCHOS_ASSERT_EQUALITY( 0, u->Scale( 0.5 ) );
      }
      // =======================================================================
      // tests involving two vectors
      Teuchos::RCP<Epetra_Vector> v = Teuchos::rcp( new Epetra_Vector( *map ) );
      v->Random();

      Teuchos::RCP<Teuchos::Time> dotTime =
          Teuchos::TimeMonitor::getNewTimer("Vector::Dot");
      {
          Teuchos::TimeMonitor tm(*dotTime);
          double dot;
          TEUCHOS_ASSERT_EQUALITY( 0, u->Dot( *v, &dot ) );
      }

      Teuchos::RCP<Teuchos::Time> multiplyTime =
          Teuchos::TimeMonitor::getNewTimer("Vector::Multiply");
      {
          Teuchos::TimeMonitor tm(*dotTime);
          TEUCHOS_ASSERT_EQUALITY( 0, u->Multiply( 1.0, *u, *v, 1.0 ) );
      }

      Teuchos::RCP<Teuchos::Time> updateTime =
          Teuchos::TimeMonitor::getNewTimer("Vector::Update");
      {
          Teuchos::TimeMonitor tm(*updateTime);
          TEUCHOS_ASSERT_EQUALITY( 0, u->Update( 1.0, *v, 1.0 ) );
      }
      // =======================================================================
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
