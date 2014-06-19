#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Epetra_LinearProblem.h>

// #include <Tpetra_Platform.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>

//#include "BelosConfigDefs.hpp"
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosMinresSolMgr.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// =============================================================================
typedef Epetra_MultiVector               MV;
typedef Epetra_Operator                  OP;
// =============================================================================
int main ( int argc, char *argv[] )
{
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm eComm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm eComm();
#endif

  Tpetra::DefaultPlatform::DefaultPlatformType &platform =
    Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > tComm = platform.getComm();

  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
  typedef Tpetra::Map<int,int,Node>                      Map;

  const Teuchos::RCP<Teuchos::FancyOStream> out =
      Teuchos::VerboseObjectBase::getDefaultOStream();

  bool success = true;
  try
  {
    // ===========================================================================
    // handle command line arguments
    Teuchos::CommandLineProcessor My_CLP;

    My_CLP.setDocString (
            "Linear solver testbed for the 1D Poisson matrix.\n"
    );

    std::string action( "matvec" );
    My_CLP.setOption("action", &action, "Which action to perform with the operator (matvec, solve_cg, solve_minres, solve_gmres)");

    std::string solver( "cg" );
//       My_CLP.setOption("solver", &solver, "Krylov subspace method (cg, minres, gmres)");

    bool verbose = true;
    My_CLP.setOption("verbose", "quiet", &verbose, "Print messages and results.");

    int frequency = 10;
    My_CLP.setOption("frequency", &frequency, "Solvers frequency for printing residuals (#iters).");

    int n = 1000;
    My_CLP.setOption("size", &n, "Size of the equation system (default: 1000).");

    // print warning for unrecognized arguments
    My_CLP.recogniseAllOptions ( true );
    My_CLP.throwExceptions ( false );

    // finally, parse the command line
    TEUCHOS_ASSERT_EQUALITY(My_CLP.parse ( argc, argv ),
                            Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL
                            );
    // =========================================================================
    // Construct Epetra matrix.
    Teuchos::RCP<Teuchos::Time> matrixConstructTime =
        Teuchos::TimeMonitor::getNewTimer("Epetra matrix construction");
    Teuchos::RCP<Epetra_CrsMatrix> epetra_A;
    {
      Teuchos::TimeMonitor tm(*matrixConstructTime);
      // Build the matrix (-1,2,-1).
      Epetra_Map map(n, 0, eComm);
      int * myGlobalElements = map.MyGlobalElements();
      epetra_A = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map, 3));
      double vals[] = {-1.0, 2.0, -1.0};
      for (int k=0; k < map.NumMyElements(); k++)
      {
        if (myGlobalElements[k] == 0)
        {
          int cols[] = {myGlobalElements[k], myGlobalElements[k]+1};
          TEUCHOS_ASSERT_EQUALITY(0, epetra_A->InsertGlobalValues(
                myGlobalElements[k],
                2,
                &vals[1],
                cols
                ));
        }
        else if (myGlobalElements[k] == n-1)
        {
          int cols[] = {myGlobalElements[k]-1, myGlobalElements[k]};
          TEUCHOS_ASSERT_EQUALITY(0, epetra_A->InsertGlobalValues(
                myGlobalElements[k],
                2,
                vals,
                cols
                ));
        }
        else
        {
          int cols[] = {myGlobalElements[k]-1, myGlobalElements[k], myGlobalElements[k]+1};
          TEUCHOS_ASSERT_EQUALITY(0, epetra_A->InsertGlobalValues(
                myGlobalElements[k],
                3,
                vals,
                cols
                ));
        }
      }
      TEUCHOS_ASSERT_EQUALITY(0, epetra_A->FillComplete(true));
    }
//       epetra_A->Print(std::cout);

    // Construct Tpetra matrix.
    Teuchos::RCP<Teuchos::Time> tpetraMatrixConstructTime =
        Teuchos::TimeMonitor::getNewTimer("Tpetra matrix construction");
    Teuchos::RCP<Tpetra::CrsMatrix<double,int> > tpetra_A;
    {
      Teuchos::TimeMonitor tm(*tpetraMatrixConstructTime);
      Teuchos::RCP<const Tpetra::Map<int> > map =
        Tpetra::createUniformContigMap<int,int>(n, tComm);
      // Get update list and number of local equations from newly created map.
      const size_t numMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const int> myGlobalElements = map->getNodeElementList();
      // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
      tpetra_A = Tpetra::createCrsMatrix<double>(map, 3);
      // Add rows one-at-a-time
      for (size_t i=0; i<numMyElements; i++)
      {
        if (myGlobalElements[i] == 0)
        {
          tpetra_A->insertGlobalValues(
              myGlobalElements[i],
              Teuchos::tuple<int>(myGlobalElements[i], myGlobalElements[i]+1),
              Teuchos::tuple<double>(2.0, -1.0)
              );
        }
        else if (myGlobalElements[i] == n-1)
        {
          tpetra_A->insertGlobalValues(
              myGlobalElements[i],
              Teuchos::tuple<int>(myGlobalElements[i]-1, myGlobalElements[i]),
              Teuchos::tuple<double>(-1.0, 2.0)
              );
        }
        else
        {
          tpetra_A->insertGlobalValues(
              myGlobalElements[i],
              Teuchos::tuple<int>(myGlobalElements[i]-1,
                                  myGlobalElements[i],
                                  myGlobalElements[i]+1
                                  ),
              Teuchos::tuple<double>(-1.0, 2.0, -1.0)
              );
        }
      }
      // Complete the fill, ask that storage be reallocated and optimized
      tpetra_A->fillComplete(Tpetra::DoOptimizeStorage);
    }
//       Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
//       tpetra_A->describe(*fos, Teuchos::VERB_EXTREME);
//       std::cout << std::endl << tpetra_A->description() << std::endl << std::endl;


    // create initial guess and right-hand side
    Teuchos::RCP<Epetra_Vector> epetra_x =
            Teuchos::rcp( new Epetra_Vector( epetra_A->OperatorDomainMap() ) );
    Teuchos::RCP<Epetra_MultiVector> epetra_b =
            Teuchos::rcp( new Epetra_Vector( epetra_A->OperatorRangeMap() ) );
    // epetra_b->Random();
    TEUCHOS_ASSERT_EQUALITY(0, epetra_b->PutScalar( 1.0 ));

    // create tpetra vectors
    Teuchos::RCP<Tpetra::Vector<double,int> > tpetra_x =
      Teuchos::rcp( new Tpetra::Vector<double,int>(tpetra_A->getDomainMap()) );
    Teuchos::RCP<Tpetra::Vector<double,int> > tpetra_b =
      Teuchos::rcp( new Tpetra::Vector<double,int>(tpetra_A->getRangeMap()) );
    tpetra_b->putScalar( 1.0 );

    if (action.compare("matvec") == 0)
    {
      TEUCHOS_ASSERT_EQUALITY(0, epetra_x->PutScalar( 1.0 ));
      Teuchos::RCP<Teuchos::Time> mvTime =
        Teuchos::TimeMonitor::getNewTimer("Epetra operator apply");
      {
        Teuchos::TimeMonitor tm(*mvTime);
        // Don't TEUCHOS_ASSERT_EQUALITY() here for speed.
        epetra_A->Apply(*epetra_x, *epetra_b);
      }

      tpetra_x->putScalar( 1.0 );
      Teuchos::RCP<Teuchos::Time> tmvTime =
        Teuchos::TimeMonitor::getNewTimer("Tpetra operator apply");
      {
        Teuchos::TimeMonitor tm(*tmvTime);
        tpetra_A->apply(*tpetra_x, *tpetra_b);
      }

      // print timing data
      Teuchos::TimeMonitor::summarize();
    }
    else
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
      Belos::LinearProblem<double,MV,OP> problem(epetra_A, epetra_x, epetra_b);
      bool set = problem.setProblem();
      TEUCHOS_TEST_FOR_EXCEPTION(
          !set,
          std::logic_error,
          "ERROR:  Belos::LinearProblem failed to set up correctly!"
          );
      // -----------------------------------------------------------------------
      // Create an iterative solver manager.
      Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver;
      if (action.compare("solve_cg") == 0)
      {
        belosList.set( "Assert Positive Definiteness", false );
        newSolver = Teuchos::rcp(
            new Belos::PseudoBlockCGSolMgr<double,MV,OP>(Teuchos::rcp(&problem,false),
            Teuchos::rcp(&belosList,false)
            ));
      }
      else if (action.compare("solve_minres") == 0)
      {
        newSolver = Teuchos::rcp(
            new Belos::MinresSolMgr<double,MV,OP>(Teuchos::rcp(&problem,false),
            Teuchos::rcp(&belosList,false)
            ));
      }
      else if (action.compare("solve_gmres") == 0)
      {
        newSolver = Teuchos::rcp(
            new Belos::PseudoBlockGmresSolMgr<double,MV,OP>(Teuchos::rcp(&problem,false),
            Teuchos::rcp(&belosList,false)
            ));
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
            true,
            "Unknown solver type \"" << solver << "\"."
            );
      }

      // Perform solve
      Teuchos::RCP<Teuchos::Time> solveTime =
        Teuchos::TimeMonitor::getNewTimer("Linear system solve");
      {
        Teuchos::TimeMonitor tm(*solveTime);
        Belos::ReturnType ret = newSolver->solve();
        success = ret==Belos::Converged;
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
