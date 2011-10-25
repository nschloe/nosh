// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_LinearProblem.h>
#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_FECrsGraph.h>

#include "Ginla_EpetraFVM_StkMeshReader.hpp"
#include "Ginla_EpetraFVM_KeoFactory.hpp"
#include "Ginla_EpetraFVM_JacobianOperator.hpp"
#include "Ginla_EpetraFVM_KeoPreconditioner.hpp"
#include "Ginla_MagneticVectorPotential_Custom.hpp"

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
//#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// =============================================================================
typedef double                           ST;
typedef Epetra_MultiVector               MV;
typedef Epetra_Operator                  OP;
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
    Anasazi::ReturnType returnCode;

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
      Teuchos::RCP<Epetra_MultiVector>    mvpValues = Teuchos::null;
      Teuchos::RCP<Epetra_Vector>         thickness = Teuchos::null;
      Teuchos::RCP<Ginla::EpetraFVM::StkMesh> mesh = Teuchos::null;

      if ( eComm->MyPID() == 0 )
          std::cout << "Reading..." << std::endl;

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

      Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp;
      double mu;
      Teuchos::RCP<Teuchos::Time> mvpConstructTime = Teuchos::TimeMonitor::getNewTimer("MVP construction");
      {
          Teuchos::TimeMonitor tm(*mvpConstructTime);
          mu = problemParameters.get<double> ( "mu" );
          mu = 2.0e-1;
          mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential::Custom ( mesh, mvpValues, mu ) );
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
          Teuchos::rcp( new Ginla::EpetraFVM::KeoFactory( mesh, thickness, mvp ) );

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
          jac = Teuchos::rcp( new Ginla::EpetraFVM::JacobianOperator( mesh, thickness, mvp, z ) );
      }

      // create preconditioner
      Teuchos::RCP<Teuchos::Time> precConstructTime = Teuchos::TimeMonitor::getNewTimer("Prec construction");
      Teuchos::RCP<Ginla::EpetraFVM::KeoPreconditioner> prec;
      if ( isPrec )
      {
          Teuchos::TimeMonitor tm(*precConstructTime);
          // create the jacobian operator
          prec = Teuchos::rcp( new Ginla::EpetraFVM::KeoPreconditioner( mesh, thickness, mvp ) );

          // actually fill it with values
          prec->rebuild();
      }


      // Create the eigensolver.
      const string which       = "LM";
      const int    nev         = 10;
      const int    blockSize   = 2;
      const int    numBlocks   = 8;
      const int    maxRestarts = 100;
      const int    maxIters    = 500;
      const double tol         = 1.0e-05;

      typedef Epetra_MultiVector MV;
      typedef Epetra_Operator OP;
      typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

      // Create an Epetra_MultiVector for an initial vector to start the solver.
      // Note:  This needs to have the same number of columns as the blocksize.
      //
      Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(jac->OperatorDomainMap(), blockSize) );
      ivec->Random();

      // Create the eigenproblem.
      //
      Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
        Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(jac, ivec) );

      // Inform the eigenproblem that the operator A is symmetric
      //
      MyProblem->setHermitian(true);

      // Set the number of eigenvalues requested
      //
      MyProblem->setNEV( nev );

      // Set the preconditioner. (May be NULL and not used.)
      MyProblem->setPrec( prec );

      // Inform the eigenproblem that you are finishing passing it information
      //
      TEUCHOS_ASSERT( MyProblem->setProblem() );

      // Create parameter list to pass into the solver manager
      //
      Teuchos::ParameterList MyPL;

      MyPL.set( "Which", which );
      MyPL.set( "Block Size", blockSize );
      MyPL.set( "Num Blocks", numBlocks );
      MyPL.set( "Maximum Restarts", maxRestarts );
      MyPL.set( "Maximum Iterations", maxIters );
      MyPL.set( "Convergence Tolerance", tol );
      MyPL.set( "Full Ortho", true );
      MyPL.set( "Use Locking", true );
      MyPL.set( "Verbosity", Anasazi::IterationDetails +
                             Anasazi::Errors +
                             Anasazi::Warnings +
                             Anasazi::FinalSummary
              );
      //
      // Create the solver manager
      Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);
//      Anasazi::BlockDavidsonSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);
//      Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);

      // Solve the problem
      //
      returnCode = MySolverMan.solve();

      // get the solution
      const Anasazi::Eigensolution<double,MV>& anasaziSolution =
              MyProblem->getSolution();

      int numVecs = anasaziSolution.numVecs;
      std::cout << "Number of computed eigenpairs: " << numVecs << std::endl;

      Teuchos::RCP<std::vector<double> > evals_r =
               Teuchos::rcp(new std::vector<double>(numVecs));
      Teuchos::RCP<std::vector<double> > evals_i =
              Teuchos::rcp(new std::vector<double>(numVecs));
      std::cout << "\n\nEigenvalues:" << std::endl;
      for (int i=0; i<numVecs; i++)
      {
          (*evals_r)[i] = anasaziSolution.Evals[i].realpart;
          (*evals_i)[i] = anasaziSolution.Evals[i].imagpart;

          std::cout << (*evals_r)[i] << " + I " << (*evals_i)[i] << std::endl;
      }
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

    return returnCode==Anasazi::Converged ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
