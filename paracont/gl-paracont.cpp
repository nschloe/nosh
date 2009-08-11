// Trilinos headers
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "NOX_Epetra_Vector.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"

// Defcont headers
#include "ContinuationManager.H"
// #include "SwiftHohenbergProblem.H"
// #include "PhaseConstraint.H"

// Main driver
int main( int argc, char **argv )
{

  try {

    // Initialise MPI
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RefCountPtr <Epetra_MpiComm> comm =
                               Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Teuchos::RefCountPtr <Epetra_SerialComm> comm =
                                          Teuchos::rcp(new Epetra_SerialComm());
#endif

    // Instantiate the continuation manager
    Teuchos::RefCountPtr<ContinuationManager> contManager =
                        Teuchos::rcp (new ContinuationManager(comm,"task.xml"));

    // Getting the task parameters
    Teuchos::RefCountPtr <Teuchos::ParameterList> taskList =
                                                     contManager->GetTaskList();

//     // Instantiate the problem
//     Teuchos::RefCountPtr <SwiftHohenbergProblem> problem =
//                          Teuchos::rcp(new SwiftHohenbergProblem(comm,taskList));
// 
//     // Set the problem in the continuation manager
//     contManager->SetLOCAProblem(problem);
// 
//     // Create the loca vector for the initial guess (needed from the constraint)
//     Teuchos::RefCountPtr <NOX::Epetra::Vector> locaInitialGuess =
//              Teuchos::rcp (new NOX::Epetra::Vector(problem->GetInitialGuess()));
// 
//     Teuchos::RefCountPtr <PhaseConstraint> phaseConstraint =
//                    Teuchos::rcp(new PhaseConstraint(problem,*locaInitialGuess));
// 
//     // Instantiate the interface
//     Teuchos::RefCountPtr <LOCA::MultiContinuation::ConstraintInterface>
//                                           interfaceConstraint = phaseConstraint;
// 
//     // Activate the constraint
//     contManager->SetLOCAConstraint(interfaceConstraint);
// 
//     // Prepare to run LOCA
//     contManager->BuildLOCAStepper();
// 
//     // Run LOCA
//     contManager->RunLOCAStepper();

  }

  catch (std::exception& e) {
    cout << e.what() << endl;
  }

  catch (const char *s) {
    cout << s << endl;
  }

  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  // Finalise MPI
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);

}