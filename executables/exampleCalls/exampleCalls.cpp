#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Vector.h>

#include "Nosh_StkMesh.hpp"
#include "Nosh_ScalarField_Constant.hpp"
#include "Nosh_MatrixBuilder_Keo.hpp"
#include "Nosh_MatrixBuilder_Laplace.hpp"
#include "Nosh_VectorField_ExplicitValues.hpp"
#include "Nosh_VectorField_ConstantCurl.hpp"
#include "Nosh_ModelEvaluator_Nls.hpp"
#include "Nosh_ModelEvaluator_Bordered.hpp"

// =============================================================================
using Teuchos::rcp;
using Teuchos::RCP;
// =============================================================================
int main(int argc, char *argv[])
{
  // Create a communicator for Epetra objects.
#ifdef HAVE_MPI
  MPI_Init( &argc, &argv );
  RCP<const Epetra_MpiComm> eComm
    = rcp<Epetra_MpiComm>(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  RCP<const Epetra_SerialComm> eComm
    = rcp<Epetra_SerialComm>(new Epetra_SerialComm());
#endif

  // Create output stream. (Handy for multicore output.)
  const RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();

  bool success = true;
  try
  {
    // ===========================================================================
    // Handle command line arguments.
    Teuchos::CommandLineProcessor myClp;

    myClp.setDocString(
      "Example usage of the model evaluator.\n"
    );

    std::string dataFile = "";
    myClp.setOption ("mesh-file", &dataFile, "Mesh file", true );

    // Print warning for unrecognized arguments and make sure to throw an
    // exception if something went wrong.
    myClp.recogniseAllOptions ( true );
    myClp.throwExceptions(true);

    // Finally, parse the command line.
    myClp.parse(argc, argv);
    // =======================================================================
    // Read the data from the file.
    const int step = 0;
    RCP<Nosh::StkMesh> mesh = rcp(new Nosh::StkMesh(*eComm, dataFile, step));

    // Cast the data into something more accessible.
    RCP<Epetra_Vector> psi = mesh->createComplexVector("psi");

    // Set the thickness field.
    RCP<Nosh::ScalarField::Virtual> thickness =
      rcp(new Nosh::ScalarField::Constant(1.0));

    // - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -
    // Some alternatives for the positive-definite operator.
    // (a) -\Delta (Laplace operator with Neumann boundary)
    //const RCP<Nosh::MatrixBuilder::Virtual> matrixBuilder =
    //  rcp(new Nosh::MatrixBuilder::Laplace(mesh, thickness));

    // (b) (-i\nabla-A)^2 (Kinetic energy of a particle in magnetic field)
    // (b1) 'A' explicitly given in file.
    const double initMu = 0.0;
    RCP<Nosh::VectorField::Virtual> mvp =
      rcp(new Nosh::VectorField::ExplicitValues(*mesh, "A", initMu));
    const RCP<Nosh::MatrixBuilder::Virtual> matrixBuilder =
      rcp(new Nosh::MatrixBuilder::Keo(mesh, thickness, mvp));

    // (b2) 'A' analytically given (here with constant curl).
    //      Optionally add a rotation axis u. This is important
    //      if continuation happens as a rotation of the vector
    //      field around an axis.
    //const RCP<DoubleVector> b = rcp(new DoubleVector(3));
    //RCP<Teuchos::SerialDenseVector<int,double> > u = Teuchos::null;
    //if ( piroParams->isSublist("Rotation vector") )
    //{
    //    u = rcp(new Teuchos::SerialDenseVector<int,double>(3));
    //    Teuchos::ParameterList & rotationVectorList =
    //        piroParams->sublist( "Rotation vector", false );
    //    (*u)[0] = rotationVectorList.get<double>("x");
    //    (*u)[1] = rotationVectorList.get<double>("y");
    //    (*u)[2] = rotationVectorList.get<double>("z");
    //}
    //RCP<Nosh::VectorField::Virtual> mvp =
    //  rcp(new Nosh::VectorField::ConstantCurl(mesh, b, u));
    //const RCP<Nosh::MatrixBuilder::Virtual> matrixBuilder =
    //  rcp(new Nosh::MatrixBuilder::Keo(mesh, thickness, mvp));
    // (b3) 'A' analytically given in a class you write yourself, derived
    //      from Nosh::MatrixBuilder::Virtual.
    // [...]
    // - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -
    // Setup the scalar potential V.
    // (a) A constant potential.
    RCP<Nosh::ScalarField::Virtual> sp =
      rcp(new Nosh::ScalarField::Constant(-1.0));
    //const double T = 0.0;
    // (b) One you built yourself by deriving from Nosh::ScalarField::Virtual.
    //RCP<Nosh::ScalarField::Virtual> sp =
      //rcp(new MyScalarField(mesh));
    // - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

    // Finally, create the model evaluator.
    // This is the most important object in the whole stack.
    const double g = 1.0;
    RCP<Nosh::ModelEvaluator::Virtual> nlsModel =
      rcp(new Nosh::ModelEvaluator::Nls(mesh, matrixBuilder, sp, g, thickness, psi));

    RCP<Nosh::ModelEvaluator::Virtual> modelEvaluator;
    const bool useBordering = false;
    if (useBordering)
    {
      // Use i*psi as bordering.
      RCP<Epetra_Vector> bordering =
        rcp(new Epetra_Vector(psi->Map()));
      for (int k=0; k<psi->Map().NumMyElements()/2; k++)
      {
        (*bordering)[2*k] = - (*psi)[2*k+1];
        (*bordering)[2*k+1] = (*psi)[2*k];
        //(*bordering)[2*k]   = 1.0;
        //(*bordering)[2*k+1] = 0.0;
      }
      //bordering->Random();
      // Initial value for the extra variable.
      double lambda = 0.0;
      modelEvaluator = rcp(new Nosh::ModelEvaluator::Bordered(nlsModel, bordering, lambda));
    }
    else
    {
      modelEvaluator = nlsModel;
    }

    // Evaluate the nonlinear Schroedinger equation F(X).
    EpetraExt::ModelEvaluator::InArgs inArgs = modelEvaluator->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs outArgs = modelEvaluator->createOutArgs();
    // Set the parameter vector.
    RCP<const Epetra_Vector> p = modelEvaluator->get_p_init(0);
    inArgs.set_p(0, p);
    // Set the X in inArgs.
    inArgs.set_x(psi);
    // Set an empty FX in outArgs. This will be filled after the call to evalModel.
    RCP<Epetra_Vector> fx = rcp(new Epetra_Vector(*modelEvaluator->get_f_map()));
    outArgs.set_f(fx);
    // Fill outArgs.
    modelEvaluator->evalModel(inArgs, outArgs);
    // fx is now filled with F(X).
    std::cout << *fx << std::endl;
    // Reset to null to make sure it's not refilled the next time evalModel is called.
    RCP<Epetra_Vector> null = Teuchos::null;
    outArgs.set_f(null);

    // Get the Jacobian.
    Teuchos::RCP<Epetra_Operator> jac = modelEvaluator->create_W();
    outArgs.set_W(jac);
    modelEvaluator->evalModel(inArgs, outArgs);
    // Now, jac contains the Jacobian operator.
    RCP<Epetra_Operator> nullOp = Teuchos::null;
    outArgs.set_W(nullOp);

    // Get the preconditioner.
    Teuchos::RCP<Epetra_Operator> prec = modelEvaluator->create_WPrec()->PrecOp;
    outArgs.set_WPrec(prec);
    modelEvaluator->evalModel(inArgs, outArgs);
    // Now, prec contains the Jacobian operator.
    outArgs.set_WPrec(nullOp);

    // Print timing data.
    Teuchos::TimeMonitor::summarize();
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