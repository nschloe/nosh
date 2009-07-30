#include <NOX.H>
#include <NOX_Epetra.H>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>

// User's application specific files 
#include "glSystem.h"
#include "stateFileReader.h"

#include <string>

using namespace std;

int main(int argc, char *argv[])
{

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
cout << "Aa"<< endl;
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
cout << "Aa"<< endl;
#else
  Epetra_SerialComm Comm;
#endif

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();

  // Check verbosity level
  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  // ---------------------------------------------------------------------------
  // set the discretization parameters
  // (possibly overwritten when reading a file
  Teuchos::ParameterList problemParameters;
  problemParameters.set("Nx",50);
  problemParameters.set("edgelength",10.0);
  problemParameters.set("H0",0.4);
  // ---------------------------------------------------------------------------


//   // ---------------------------------------------------------------------------
//   std::vector<double_complex> psi;
//   StateFileReader fileReader;
//   std::string file = "data/solution.vtk";
//   std::string format = "legacyVTK";
//   fileReader.readFile( file, format, &psi, &problemParameters );
//   // ---------------------------------------------------------------------------


  // create the gl problem
  GinzburgLandau glProblem = GinzburgLandau( problemParameters.get<int>("Nx"),
                                             problemParameters.get<double>("edgelength"),
                                             problemParameters.get<double>("H0")
                                           );

  // Create the interface between NOX and the application
  // This object is derived from NOX::Epetra::Interface
  Teuchos::RCP<GlSystem> glsystem = Teuchos::rcp(new GlSystem( glProblem,
                                                               Comm )
                                                );

  // Get the vector from the Problem
  Teuchos::RCP<Epetra_Vector> soln = glsystem->getSolution();
  Teuchos::RCP<NOX::Epetra::Vector> noxSoln =
    Teuchos::rcp(new NOX::Epetra::Vector(soln,
                                         NOX::Epetra::Vector::CreateView));

  // Set the initial guess 
//  soln->PutScalar(0.5);

  // Begin Nonlinear Solver ************************************
  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID); 
  printParams.set("Output Precision", 16);
  printParams.set("Output Processor", 0);
  if (verbose)
    printParams.set("Output Information", 
                             NOX::Utils::OuterIteration +
                             NOX::Utils::OuterIterationStatusTest +
                             NOX::Utils::InnerIteration +
                             NOX::Utils::LinearSolverDetails +
                             NOX::Utils::Parameters +
                             NOX::Utils::Details +
                             NOX::Utils::Warning +
                             NOX::Utils::Debug +
                             NOX::Utils::TestDetails +
                             NOX::Utils::Error );
  else
    printParams.set("Output Information", NOX::Utils::Error +
                             NOX::Utils::TestDetails);

  // Create a print class for controlling output below
  NOX::Utils printing(printParams);

  // Sublist for line search

  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
//  lsParams.set("Amesos Solver", "Amesos_Superlu");

  // lsParams.set("Aztec Solver", "BiCGStab");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Output Frequency", 32);
  lsParams.set("Output Solver Details", true);
  lsParams.set("Max Iterations", 2000);
  lsParams.set("Tolerance", 1e-4);


  // ====================== BEGIN IFPACK PRECONDITIONER ========================
  lsParams.set("Preconditioner","New Ifpack");
  lsParams.set("Ifpack Preconditioner","Amesos");

  lsParams.set("Use Preconditioner as Solver",true);
  lsParams.set("Max Age Of Prec", 1);

  Teuchos::ParameterList& IFPACKparams = lsParams.sublist("Ifpack");
      IFPACKparams.set("amesos: solver type","Amesos_Superlu");
      IFPACKparams.set("fact: ilut level-of-fill",1.0);
      IFPACKparams.set("fact: drop tolerance", 1e-1);
      IFPACKparams.set("schwarz: combine mode", "Zero");
      IFPACKparams.set("schwarz: compute condest", true);
  // ====================== END IFPACK PRECONDITIONER ==========================

  // Let's force all status tests to do a full check
  nlParams.sublist("Solver Options").set("Status Test Check Type", "Complete");

  // Create all possible Epetra_Operators.
  Teuchos::RCP<Epetra_RowMatrix> Analytic = glsystem->getJacobian();

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = glsystem;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glsystem;

  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
           Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO( printParams,
                                                              lsParams,
                                                              iReq,
                                                              iJac,
                                                              Analytic,
                                                              *soln        ) );

  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::Epetra::Vector::CreateView);
  Teuchos::RCP<NOX::Epetra::Group> grpPtr =
    Teuchos::rcp(new NOX::Epetra::Group(printParams,
                                        iReq,
                                        initialGuess,
                                        linSys));

  NOX::Epetra::Group& grp = *grpPtr;
  // ------------------------------------------------------------------------


  // ------------------------------------------------------------------------
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-12));
  Teuchos::RCP<NOX::StatusTest::NormF> relresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-12));
  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-12));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-5));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(relresid);
  converged->addStatusTest(wrms);
  converged->addStatusTest(update);

  int maxNonlinearIterations = 20;
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNonlinearIterations));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  // this test is useful if we start in a solution
  Teuchos::RCP<NOX::StatusTest::NormF> absresexact =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-20));
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  combo->addStatusTest(fv);
  combo->addStatusTest(absresexact);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);
  // ------------------------------------------------------------------------

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
  NOX::StatusTest::StatusType        solvStatus = solver->solve();
  // End Nonlinear Solver **************************************

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = 
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& finalSolution = 
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).
    getEpetraVector();

//   // get the energy
//   passVector sol;
//   sol.values = finalSolution.Values();
//   sol.length = soln->Map().NumMyElements();
//   get_energy_libgl(&sol,&energy);
//   cout << "energy:"<< energy << endl;

  // Output the parameter list
  if (verbose) {
    if (printing.isPrintType(NOX::Utils::Parameters)) {
      printing.out() << endl << "Final Parameters" << endl
           << "****************" << endl;
      solver->getList().print(printing.out());
      printing.out() << endl;

    }
  }

  // Print solution to a file
  char file_name[25];
  FILE *ifp;
  int NumMyElements = soln->Map().NumMyElements();
  (void) sprintf(file_name, "output.%d",MyPID);
  ifp = fopen(file_name, "w");
  for ( int i=0; i<NumMyElements; i++ )
    fprintf( ifp, "%d  %E\n", soln->Map().MinMyGID()+i, finalSolution[i] );
  fclose(ifp);

  // ---------------------------------------------------------------------------
  // print the solution to a file
  std::string  fileFormat = "XDMF";
  std::string  fileName   = "data/solution.xmf";
  glsystem->solutionToFile( finalSolution,
                            problemParameters,
                            fileFormat,
                            fileName  );
  // ---------------------------------------------------------------------------


  // gather full solution to one processor, put it to a file
  nlParams.print(cout,1,true,true);

  // Tests
  int status = 0; // Converged

  // 1. Convergence
  if (solvStatus != NOX::StatusTest::Converged) {
      status = 1;
      if (printing.isPrintType(NOX::Utils::Error))
        printing.out() << "Nonlinear solver failed to converge!" << endl;
  }
#ifndef HAVE_MPI
  // 2. Linear solve iterations (53) - SERIAL TEST ONLY!
  //    The number of linear iterations changes with # of procs.
  if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").get("Total Number of Linear Iterations",0) != 53) {
    status = 2;
  }
#endif
  // 3. Nonlinear solve iterations (10)
  if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Output").get("Nonlinear Iterations", 0) == maxNonlinearIterations)
    status = 3;

  // 

//   // 4. Test the pre/post iterate options
//   {
//   UserPrePostOperator* ppoPtr = dynamic_cast<UserPrePostOperator*>(ppo.get());
//   if (ppoPtr->getNumRunPreIterate() != 10)
//     status = 4;
//   if (ppoPtr->getNumRunPostIterate() != 10)
//     status = 4;
//   if (ppoPtr->getNumRunPreSolve() != 1)
//     status = 4;
//   if (ppoPtr->getNumRunPostSolve() != 1)
//     status = 4;
//   }
  // Summarize test results 
  if (status == 0)
    printing.out() << "Test passed!" << endl;
  else 
    printing.out() << "Test failed!" << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  printing.out() << "Status = " << status << endl;

cout << "Final parameters:" << endl
     << "===================" << endl
     << solver->getList();

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
