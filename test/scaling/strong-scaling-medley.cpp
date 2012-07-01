// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

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
      // Create map.
      // Do strong scaling tests, so keep numGlobalElements independent of
      // the number of processes.
      int numGlobalElements = 5e7;
      int indexBase = 0;
      Teuchos::RCP<Epetra_Map> map =
          Teuchos::rcp( new Epetra_Map ( numGlobalElements, indexBase, *eComm ) );
      // =======================================================================
//       // Create map with overlay.
//       int numMyOverlapNodes = 3;
// 
//       // Get an approximation of my nodes.
//       int numMyElements = numGlobalElements / eComm->NumProc();
//       int startIndex = eComm->MyPID() * numMyElements;
//       // Calculate the resulting number of total nodes.
//       int numTotalNodes = numMyElements * eComm->NumProc();
//       // Add one node to the first numGlobalElements-numTotalNodes processes.
//       if ( eComm->MyPID() < numGlobalElements - numTotalNodes )
//       {
//           numMyElements++;
//           startIndex += eComm->MyPID();
//       }
//       else
//       {
//           startIndex += numGlobalElements - numTotalNodes;
//       }
// 
//       Teuchos::Array<int> indices( numMyElements );
//       for ( int k = 0;  k<numMyElements; k++ )
//           indices[k] = startIndex + k;
// 
//       std::cout << numGlobalElements << std::endl;
//       std::cout << numMyElements << std::endl;
// 
//       Teuchos::RCP<Epetra_Map> overlapMap =
//           Teuchos::rcp( new Epetra_Map ( numGlobalElements, numMyElements, indices.getRawPtr(), indexBase, *eComm ) );
// 
//       overlapMap->Print( std::cout );
// 
//       throw 1;
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
          Teuchos::TimeMonitor tm(*multiplyTime);
          TEUCHOS_ASSERT_EQUALITY( 0, u->Multiply( 1.0, *u, *v, 1.0 ) );
      }

      Teuchos::RCP<Teuchos::Time> updateTime =
          Teuchos::TimeMonitor::getNewTimer("Vector::Update");
      {
          Teuchos::TimeMonitor tm(*updateTime);
          TEUCHOS_ASSERT_EQUALITY( 0, u->Update( 1.0, *v, 1.0 ) );
      }
      // =======================================================================
      // matrix-vector tests

      // diagonal test matrix
      Teuchos::RCP<Epetra_CrsMatrix> D =
          Teuchos::rcp( new Epetra_CrsMatrix( Copy, *map, 1 ) );
      for ( int k=0; k < map->NumMyElements(); k++ )
      {
          int col = map->GID(k);
          double val = 1.0 / (col+1);
//          TEUCHOS_ASSERT_EQUALITY( 0, D->InsertMyValues( k, 1, &val, &col ) );
          TEUCHOS_ASSERT_EQUALITY( 0, D->InsertGlobalValues( col, 1, &val, &col ) );
      }
      TEUCHOS_ASSERT_EQUALITY( 0, D->FillComplete() );

      // tridiagonal test matrix
      Teuchos::RCP<Epetra_CrsMatrix> T =
          Teuchos::rcp( new Epetra_CrsMatrix( Copy, *map, 3 ) );
      for ( int k=0; k < map->NumMyElements(); k++ )
      {
          int row = map->GID(k);
          if ( row > 0 )
          {
              int col = row-1;
              double val = -1.0;
//              TEUCHOS_ASSERT_EQUALITY( 0, T->InsertMyValues( k, 1, &val, &col ) );
              TEUCHOS_ASSERT_EQUALITY( 0, T->InsertGlobalValues( row, 1, &val, &col ) );
          }
          {
              int col = row;
              double val = 2.0;
//              TEUCHOS_ASSERT_EQUALITY( 0, T->InsertMyValues( k, 1, &val, &col ) );
              TEUCHOS_ASSERT_EQUALITY( 0, T->InsertGlobalValues( row, 1, &val, &col ) );
          }
          if ( row < numGlobalElements-1 )
          {
              int col = row+1;
              double val = -1.0;
//              TEUCHOS_ASSERT_EQUALITY( 0, T->InsertMyValues( k, 1, &val, &col ) );
              TEUCHOS_ASSERT_EQUALITY( 0, T->InsertGlobalValues( row, 1, &val, &col ) );
          }

      }
      TEUCHOS_ASSERT_EQUALITY( 0, T->FillComplete() );

      // start timings
      Teuchos::RCP<Teuchos::Time> mNorm1Time =
          Teuchos::TimeMonitor::getNewTimer("CrsMatrix::Norm1");
      {
          Teuchos::TimeMonitor tm(*mNorm1Time);
          double dNorm1 = D->NormOne();
          double tNorm1 = T->NormOne();
      }

      Teuchos::RCP<Teuchos::Time> mNormInfTime =
          Teuchos::TimeMonitor::getNewTimer("CrsMatrix::NormInf");
      {
          Teuchos::TimeMonitor tm(*mNormInfTime);
          double dNormInf = D->NormInf();
          double tNormInf = T->NormInf();
      }

      Teuchos::RCP<Teuchos::Time> mNormFrobTime =
          Teuchos::TimeMonitor::getNewTimer("CrsMatrix::NormFrobenius");
      {
          Teuchos::TimeMonitor tm(*mNormFrobTime);
          double dNormFrob = D->NormFrobenius();
          double tNormFrob = T->NormFrobenius();
      }

      Teuchos::RCP<Teuchos::Time> mScaleTime =
          Teuchos::TimeMonitor::getNewTimer("CrsMatrix::Scale");
      {
          Teuchos::TimeMonitor tm(*mScaleTime);
          TEUCHOS_ASSERT_EQUALITY( 0, D->Scale( 2.0 ) );
          TEUCHOS_ASSERT_EQUALITY( 0, T->Scale( 2.0 ) );
      }

      Teuchos::RCP<Teuchos::Time> leftScaleTime =
          Teuchos::TimeMonitor::getNewTimer("CrsMatrix::LeftScale");
      {
          Teuchos::TimeMonitor tm(*leftScaleTime);
          TEUCHOS_ASSERT_EQUALITY( 0, D->LeftScale( *v ) );
          TEUCHOS_ASSERT_EQUALITY( 0, T->LeftScale( *v ) );
      }

      Teuchos::RCP<Teuchos::Time> rightScaleTime =
          Teuchos::TimeMonitor::getNewTimer("CrsMatrix::RightScale");
      {
          Teuchos::TimeMonitor tm(*rightScaleTime);
          TEUCHOS_ASSERT_EQUALITY( 0, D->RightScale( *v ) );
          TEUCHOS_ASSERT_EQUALITY( 0, T->RightScale( *v ) );
      }

      Teuchos::RCP<Teuchos::Time> applyTime =
          Teuchos::TimeMonitor::getNewTimer("CrsMatrix::Apply");
      {
          Teuchos::TimeMonitor tm(*applyTime);
          TEUCHOS_ASSERT_EQUALITY( 0, D->Apply( *u, *v ) );
          TEUCHOS_ASSERT_EQUALITY( 0, T->Apply( *u, *v ) );
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
