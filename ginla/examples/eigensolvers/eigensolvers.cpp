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
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Epetra_CrsGraph.h>

#include "Ginla_StkMeshReader.hpp"
#include "Ginla_KeoFactory.hpp"
#include "Ginla_JacobianOperator.hpp"
#include "Ginla_KeoRegularized.hpp"
#include "Ginla_MagneticVectorPotential.hpp"

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

      std::string inputFilePath( "" );
      My_CLP.setOption ( "input", &inputFilePath, "Input state file", true );

      bool verbose = true;
      My_CLP.setOption("verbose","quiet",&verbose,"Print messages and results.");

      bool isPrec = true;
      My_CLP.setOption("prec","noprec",&isPrec,"Use a preconditioner.");

      int frequency = 10;
      My_CLP.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");

      std::string method = "lobpcg";
      My_CLP.setOption("method",&method,"Method for solving the eigenproblem. (*lobpcg, krylovschur, davidson)");


      // print warning for unrecognized arguments
      My_CLP.recogniseAllOptions ( true );

      // finally, parse the command line
      TEUCHOS_ASSERT_EQUALITY( My_CLP.parse ( argc, argv ),
                               Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL
                             );
      // =========================================================================
      // Read the data from the file.
      Teuchos::ParameterList data;
      Ginla::StkMeshRead( *eComm, inputFilePath, data );

      // Cast the data into something more accessible.
      Teuchos::RCP<Ginla::StkMesh>     & mesh = data.get<Teuchos::RCP<Ginla::StkMesh> >( "mesh" );
      Teuchos::RCP<Epetra_Vector>      & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );
      Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::ParameterList           & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

      Teuchos::RCP<Ginla::MagneticVectorPotential> mvp;
      double mu;
      Teuchos::RCP<Teuchos::Time> mvpConstructTime = Teuchos::TimeMonitor::getNewTimer("MVP construction");
      {
          Teuchos::TimeMonitor tm(*mvpConstructTime);
          mu = problemParameters.get<double> ( "mu" );
          mu = 2.0e-1;
          mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential ( mesh, mvpValues, mu ) );
      }

      Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
          Teuchos::rcp( new LOCA::ParameterVector() );
      mvpParameters->addParameter( "mu", mu );

      Teuchos::RCP<Ginla::KeoFactory> keoFactory =
          Teuchos::rcp( new Ginla::KeoFactory( mesh, thickness, mvp ) );

      // create Jacobian
      Teuchos::RCP<Teuchos::Time> jacobianConstructTime = Teuchos::TimeMonitor::getNewTimer("Jacobian construction");
      Teuchos::RCP<Ginla::JacobianOperator> jac;
      {
          Teuchos::TimeMonitor tm(*jacobianConstructTime);
          // create the jacobian operator
          jac = Teuchos::rcp( new Ginla::JacobianOperator( mesh, thickness, keoFactory, z ) );
      }

      // create preconditioner
      Teuchos::RCP<Teuchos::Time> precConstructTime = Teuchos::TimeMonitor::getNewTimer("Prec construction");
      Teuchos::RCP<Ginla::KeoRegularized> keoReg;
      if ( isPrec )
      {
          Teuchos::TimeMonitor tm(*precConstructTime);
          // create the jacobian operator
          keoReg = Teuchos::rcp( new Ginla::KeoRegularized( keoFactory ) );

          // actually fill it with values
          keoReg->rebuild();
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
      MyProblem->setPrec( keoReg );

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

      // Create the solver manager and solve the problem.
      //
      Anasazi::ReturnType returnCode;
      if ( method.compare("lobpcg") == 0 )
      {
          Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);
          returnCode = MySolverMan.solve();
      }
      else if ( method.compare("davidson") == 0 )
      {
          Anasazi::BlockDavidsonSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);
          returnCode = MySolverMan.solve();
      }
      else if ( method.compare("krylovschur") == 0 )
      {
          Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);
          returnCode = MySolverMan.solve();
      }
      else
          TEST_FOR_EXCEPT_MSG( true, "Invalid eigensolver method \"" << method << "\"." );

      // Check for success.
      success = returnCode==Anasazi::Converged;

      // get the solution
      const Anasazi::Eigensolution<double,MV>& anasaziSolution =
              MyProblem->getSolution();

      int numVecs = anasaziSolution.numVecs;
      *out << "Number of computed eigenpairs: " << numVecs << std::endl;

      Teuchos::ArrayRCP<double> evals_r( numVecs );
      Teuchos::ArrayRCP<double> evals_i( numVecs );
      *out << "\n\nEigenvalues:" << std::endl;
      for (int i=0; i<numVecs; i++)
      {
          evals_r[i] = anasaziSolution.Evals[i].realpart;
          evals_i[i] = anasaziSolution.Evals[i].imagpart;

          *out << evals_r[i] << " + I " << evals_i[i] << std::endl;
      }
      // -----------------------------------------------------------------------
      // print timing data
      Teuchos::TimeMonitor::summarize();
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

#ifdef HAVE_MPI
      MPI_Finalize();
#endif

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
