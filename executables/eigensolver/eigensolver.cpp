#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_LinearProblem.h>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Epetra_CrsGraph.h>

#include <AnasaziConfigDefs.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziLOBPCGSolMgr.hpp>
#include <AnasaziBlockDavidsonSolMgr.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
//#include <AnasaziBasicOutputManager.hpp>
#include <AnasaziEpetraAdapter.hpp>

#include "Nosh_MatrixBuilder_Keo.hpp"
#include "Nosh_ScalarField_Constant.hpp"
#include "Nosh_VectorField_ExplicitValues.hpp"
#include "Nosh_ModelEvaluator_Nls.hpp"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// =============================================================================
typedef double             ST;
typedef Epetra_MultiVector MV;
typedef Epetra_Operator    OP;
typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;
// =============================================================================
using Teuchos::rcp;
using Teuchos::RCP;
// =============================================================================
int main ( int argc, char *argv[] )
{
  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init( &argc, &argv );
#endif

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  RCP<Epetra_MpiComm> eComm =
    rcp<Epetra_MpiComm>(new Epetra_MpiComm (MPI_COMM_WORLD));
#else
  RCP<Epetra_SerialComm> eComm =
    rcp<Epetra_SerialComm>(new Epetra_SerialComm());
#endif

  const RCP<Teuchos::FancyOStream> out =
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

    std::string inputFilePath("");
    My_CLP.setOption("input", &inputFilePath, "input state file", true);

    std::string mvpFilePath("");
    My_CLP.setOption("mvpfile",&mvpFilePath,"file containing magnetic vector potential (default: same as input file)");

    int step = 0;
    My_CLP.setOption("step",&step,"step to read");

    bool verbose = true;
    My_CLP.setOption("verbose","quiet",&verbose,"print messages and results");

    bool isPrec = true;
    My_CLP.setOption("prec","noprec",&isPrec,"use a preconditioner");

    int frequency = 10;
    My_CLP.setOption("frequency",&frequency,"solver frequency for printing eigenresiduals (#iters)");

    double mu = 0.0;
    My_CLP.setOption("mu", &mu, "parameter value mu");

    std::string method = "lobpcg";
    My_CLP.setOption("method",&method,"method for solving the eigenproblem {lobpcg, krylovschur, davidson}");

    int numEv = 10;
    My_CLP.setOption("numev",&numEv,"number of eigenvalues to compute");

    int numBlocks = 10;
    My_CLP.setOption("numblocks",&numBlocks,"number of blocks");

    int blockSize = 2;
    My_CLP.setOption("blocksize",&blockSize,"block size");

    int maxRestarts = 10;
    My_CLP.setOption("maxrestarts",&maxRestarts,"maximum number of restarts");

    int maxIter = 100;
    My_CLP.setOption("maxiter",&maxIter,"maximum number of iterations");

    double tol = 1.0e-5;
    My_CLP.setOption("tolerance",&tol,"tolerance");

    // print warning for unrecognized arguments
    My_CLP.recogniseAllOptions( true );

    // finally, parse the command line
    TEUCHOS_ASSERT_EQUALITY(My_CLP.parse(argc, argv),
                            Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL);
    // =========================================================================
    // Read the data from the file.
    RCP<Nosh::StkMesh> mesh = rcp(new Nosh::StkMesh(*eComm, inputFilePath, step));

    const double time = mesh->getTime();
    mu = time;

    // Cast the data into something more accessible.
    RCP<Epetra_Vector> psi = mesh->createComplexVector("psi");

    // Create MVP (possibly from another mesh file).
    RCP<Nosh::VectorField::Virtual> mvp;
    RCP<Teuchos::Time> mvpConstructTime =
      Teuchos::TimeMonitor::getNewTimer("MVP construction");
    {
      Teuchos::TimeMonitor tm(*mvpConstructTime);
      if (mvpFilePath.empty())
        mvp = rcp(new Nosh::VectorField::ExplicitValues(*mesh, "A", mu));
      else
      {
        Nosh::StkMesh mesh2(*eComm, mvpFilePath, 0);
        mvp = rcp(new Nosh::VectorField::ExplicitValues(mesh2, "A", mu));
      }
    }

    // Construct thickness.
    RCP<Nosh::ScalarField::Virtual> thickness =
      rcp(new Nosh::ScalarField::Constant(1.0));

    // Create matrix builder.
    const RCP<Nosh::MatrixBuilder::Virtual> matrixBuilder =
      rcp(new Nosh::MatrixBuilder::Keo(mesh, thickness, mvp));

    // Construct scalar potential.
    RCP<Nosh::ScalarField::Virtual> sp =
      rcp(new Nosh::ScalarField::Constant(-1.0));

    // Finally, create the model evaluator.
    // This is the most important object in the whole stack.
    const double g = 1.0;
    RCP<Nosh::ModelEvaluator::Virtual> modelEvaluator =
      rcp(new Nosh::ModelEvaluator::Nls(mesh, matrixBuilder, sp, g, thickness, psi));

    // Set the input arguments.
    EpetraExt::ModelEvaluator::InArgs inArgs = modelEvaluator->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs outArgs = modelEvaluator->createOutArgs();
    RCP<const Epetra_Vector> p = modelEvaluator->get_p_init(0);
    inArgs.set_p(0, p);
    inArgs.set_x(psi);

    // Get the preconditioner.
    modelEvaluator->evalModel(inArgs, outArgs);

    // Create and fill Jacobian.
    RCP<Teuchos::Time> jacobianConstructTime =
      Teuchos::TimeMonitor::getNewTimer("Operator construction");
    RCP<Epetra_Operator> jac;
    RCP<Epetra_Operator> prec;
    // Fill outArgs.
    modelEvaluator->evalModel(inArgs, outArgs);
    {
      Teuchos::TimeMonitor tm(*jacobianConstructTime);
      // Get the Jacobian.
      jac = modelEvaluator->create_W();
      outArgs.set_W(jac);
      if (isPrec)
      {
        prec = modelEvaluator->create_WPrec()->PrecOp;
        outArgs.set_WPrec(prec);
      }
      modelEvaluator->evalModel(inArgs, outArgs);
      // Now, jac and prec contain the operators.
      // Reset outArgs.
      outArgs.set_W(Teuchos::null);
      outArgs.set_WPrec(Teuchos::null);
    }

    bool checkFx = true;
    if (checkFx)
    {
      RCP<Epetra_Vector> fx =
        rcp(new Epetra_Vector(*modelEvaluator->get_f_map()));
      outArgs.set_f(fx);
      modelEvaluator->evalModel(inArgs, outArgs);
      // Check ||F(x)||.
      double r;
      fx->Norm2(&r);
      *out << "||F(x)|| = " << r << std::endl;
    }

    // Create the eigensolver.
    // Create an Epetra_MultiVector for an initial vector to start the solver.
    // Note:  This needs to have the same number of columns as the blocksize.
    RCP<Epetra_MultiVector> ivec =
      rcp(new Epetra_MultiVector(jac->OperatorDomainMap(), blockSize));
    ivec->Random();

    // Create the eigenproblem.
    RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
      rcp(new Anasazi::BasicEigenproblem<double, MV, OP>(jac, ivec));

    // Inform the eigenproblem that the operator A is symmetric
    MyProblem->setHermitian(true);

    // Set the number of eigenvalues requested
    MyProblem->setNEV(numEv);

    // Set the preconditioner. (May be NULL and not used.)
    MyProblem->setPrec( prec );

    // Inform the eigenproblem that you are finishing passing it information
    TEUCHOS_ASSERT( MyProblem->setProblem() );

    // Create parameter list to pass into the solver manager
    Teuchos::ParameterList MyPL;

    MyPL.set("Which", "LM");
    MyPL.set("Block Size", blockSize);
    MyPL.set("Num Blocks", numBlocks);
    MyPL.set("Maximum Restarts", maxRestarts);
    MyPL.set("Maximum Iterations", maxIter);
    MyPL.set("Convergence Tolerance", tol);
    MyPL.set("Full Ortho", true);
    MyPL.set("Use Locking", true);
    MyPL.set("Verbosity", Anasazi::IterationDetails +
                          Anasazi::Errors +
                          Anasazi::Warnings +
                          Anasazi::StatusTestDetails +
                          Anasazi::Debug +
                          Anasazi::FinalSummary
                          );

    // Create the solver manager and solve the problem.
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
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                                  "Invalid eigensolver method \"" << method << "\"." );

    // Check for success.
    success = returnCode==Anasazi::Converged;

    // get the solution
    const Anasazi::Eigensolution<double,MV>& anasaziSolution =
      MyProblem->getSolution();

    const int numVecs = anasaziSolution.numVecs;
    *out << "Number of computed eigenpairs: " << numVecs << std::endl;

    Teuchos::ArrayRCP<double> evals_r( numVecs );
    Teuchos::ArrayRCP<double> evals_i( numVecs );
    if (numVecs > 0)
    {
      *out << "\nEigenvalues:" << std::endl;
      for (int i=0; i<numVecs; i++)
      {
        evals_r[i] = anasaziSolution.Evals[i].realpart;
        evals_i[i] = anasaziSolution.Evals[i].imagpart;

        *out << evals_r[i] << " + I " << evals_i[i] << std::endl;
      }
      // Check residuals
    }
    // -----------------------------------------------------------------------
    // print timing data
    //Teuchos::TimeMonitor::summarize();
  }
  catch (Teuchos::CommandLineProcessor::HelpPrinted)
  {}
  catch (Teuchos::CommandLineProcessor::ParseError)
  {}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
