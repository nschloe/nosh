#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Vector.h>

#include <Piro_Epetra_NOXSolver.hpp>
#include <Piro_Epetra_LOCASolver.hpp>

#include "Nosh_StkMesh.hpp"
#include "Nosh_ScalarField_Constant.hpp"
#include "Nosh_MatrixBuilder_Keo.hpp"
#include "Nosh_MatrixBuilder_Laplace.hpp"
#include "Nosh_VectorField_ExplicitValues.hpp"
#include "Nosh_VectorField_ConstantCurl.hpp"
#include "Nosh_ModelEvaluator_Nls.hpp"
#include "Nosh_ModelEvaluator_Bordered.hpp"
#include "Nosh_Observer.hpp"
#include "Nosh_SaveEigenData.hpp"
#include "Nosh_CsvWriter.hpp"

#include "MyScalarField.hpp"

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

  // Wrap the whole code in a big try-catch-statement.
  bool success = true;
  try
  {
    // ===========================================================================
    // Handle command line arguments.
    Teuchos::CommandLineProcessor myClp;

    myClp.setDocString(
      "This program solves nonlinear Schrödinger equations with a Piro interface.\n"
    );

    std::string xmlInputPath = "";
    myClp.setOption ("xml-input-file", &xmlInputPath,
                      "XML file containing the parameter list", true );

    // Print warning for unrecognized arguments and make sure to throw an
    // exception if something went wrong.
    myClp.recogniseAllOptions ( true );
    myClp.throwExceptions(true);

    // Finally, parse the command line.
    myClp.parse(argc, argv);

    // Retrieve Piro parameter list from given file.
    RCP<Teuchos::ParameterList> piroParams =
        rcp(new Teuchos::ParameterList());
    Teuchos::updateParametersFromXmlFile(xmlInputPath,
                                         piroParams.ptr());
    // =======================================================================
    // Extract the location of input and output files.
    const Teuchos::ParameterList outputList = piroParams->sublist("Output", true);

    // Set default directory to be the directory of the XML file itself
    const std::string xmlDirectory =
      xmlInputPath.substr(0, xmlInputPath.find_last_of( "/" ) + 1);

    const std::string & outputDirectory = xmlDirectory;

    const std::string contFilePath = xmlDirectory + "/"
                                   + outputList.get<std::string>( "Continuation data file name" );

    Teuchos::ParameterList & inputDataList =
      piroParams->sublist ( "Input", true );
    const std::string dataFile = xmlDirectory + "/"
                               + inputDataList.get<std::string>( "File" );
    const int step = inputDataList.get<int>("Initial Psi Step");

    const bool useBordering = piroParams->get<bool>("Bordering");
    // =======================================================================
    // Get the initial parameter values.
    Teuchos::ParameterList initialParameterValues =
      piroParams->sublist("Initial parameter values", true);

    // Read the data from the file.
    RCP<Nosh::StkMesh> mesh = rcp(new Nosh::StkMesh(*eComm, dataFile, 0));

    // Cast the data into something more accessible.
    RCP<Epetra_Vector> psi = mesh->createComplexVector("psi");

    // Set the output directory for later plotting with this.
    mesh->openOutputChannel(outputDirectory, "solution");

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
    const double initMu = initialParameterValues.get<double>("mu");
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
    //const double T = initialParameterValues.get<double>("T");
    // (b) One you built yourself by deriving from Nosh::ScalarField::Virtual.
    //RCP<Nosh::ScalarField::Virtual> sp =
      //rcp(new MyScalarField(mesh));
    // - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

    // Finally, create the model evaluator.
    // This is the most important object in the whole stack.
    const double g = initialParameterValues.get<double>("g");
    RCP<Nosh::ModelEvaluator::Virtual> nlsModel =
      rcp(new Nosh::ModelEvaluator::Nls(mesh, matrixBuilder, sp, g, thickness, psi));

    RCP<Nosh::ModelEvaluator::Virtual> modelEvaluator;
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

    // Build the Piro model evaluator. It's used to hook up with
    // several different backends (NOX, LOCA, Rhythmos,...).
    RCP<EpetraExt::ModelEvaluator> piro;

    // Declare the eigensaver; it will be used only for LOCA solvers, though.
    RCP<Nosh::SaveEigenData> glEigenSaver;

    // Switch by solver type.
    std::string & solver = piroParams->get<std::string>("Piro Solver");
    // ----------------------------------------------------------------------
    if (solver == "NOX")
    {
      RCP<Nosh::Observer> observer =
        rcp(new Nosh::Observer(modelEvaluator,
                               contFilePath,
                               Nosh::Observer::OBSERVER_TYPE_NEWTON));

      piro = rcp(new Piro::Epetra::NOXSolver(piroParams,
                                             modelEvaluator,
                                             observer));
    }
    // ----------------------------------------------------------------------
    else if (solver == "LOCA")
    {
      RCP<Nosh::Observer> observer =
        rcp(new Nosh::Observer(modelEvaluator,
                               contFilePath,
                               Nosh::Observer::OBSERVER_TYPE_CONTINUATION,
                               piroParams->sublist("LOCA")
                                          .sublist("Stepper")
                                          .get<std::string>("Continuation Parameter")
                               ));

      // Setup eigen saver.
#ifdef HAVE_LOCA_ANASAZI
      bool computeEigenvalues = piroParams->sublist( "LOCA" )
                                           .sublist( "Stepper" )
                                           .get<bool>("Compute Eigenvalues");
      if (computeEigenvalues)
      {
        Teuchos::ParameterList & eigenList = piroParams->sublist("LOCA")
                                                        .sublist("Stepper")
                                                        .sublist("Eigensolver");
        std::string eigenvaluesFilePath = xmlDirectory
                                        + "/"
                                        + outputList.get<std::string> ( "Eigenvalues file name" );

        glEigenSaver =
          RCP<Nosh::SaveEigenData>(new Nosh::SaveEigenData(eigenList,
                                                           modelEvaluator,
                                                           eigenvaluesFilePath));

        RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy =
          glEigenSaver;
        eigenList.set("Save Eigen Data Method",
                      "User-Defined");
        eigenList.set("User-Defined Save Eigen Data Name",
                      "glSaveEigenDataStrategy");
        eigenList.set("glSaveEigenDataStrategy",
                      glSaveEigenDataStrategy);
      }
#endif
      // Get the solver.
      RCP<Piro::Epetra::LOCASolver> piroLOCASolver =
        rcp(new Piro::Epetra::LOCASolver(piroParams, modelEvaluator, observer));

      // Get stepper and inject it into the eigensaver.
      RCP<LOCA::Stepper> stepper = piroLOCASolver->getLOCAStepperNonConst();
#ifdef HAVE_LOCA_ANASAZI
      if (computeEigenvalues)
        glEigenSaver->setLocaStepper( stepper );
#endif
      piro = piroLOCASolver;
    }
    // ----------------------------------------------------------------------
    else if ( solver == "Turning Point" )
    {
      RCP<Nosh::Observer> observer = Teuchos::null;

      Teuchos::ParameterList & bifList =
        piroParams->sublist("LOCA").sublist("Bifurcation");

      // Fetch the (approximate) null state.
      RCP<Epetra_Vector> nullstateZ = mesh->createVector("null");

      // Set the length normalization vector to be the initial null vector.
      TEUCHOS_ASSERT( !nullstateZ.is_null() );
      RCP<NOX::Abstract::Vector> lengthNormVec =
        rcp(new NOX::Epetra::Vector(*nullstateZ));
      //lengthNormVec->init(1.0);
      bifList.set("Length Normalization Vector", lengthNormVec);

      // Set the initial null vector.
      RCP<NOX::Abstract::Vector> initialNullAbstractVec =
        rcp(new NOX::Epetra::Vector(*nullstateZ));
      // initialNullAbstractVec->init(1.0);
      bifList.set("Initial Null Vector", initialNullAbstractVec);

      piro = rcp(new Piro::Epetra::LOCASolver(piroParams,
                                              modelEvaluator,
                                              observer));
    }
    // ----------------------------------------------------------------------
    else
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                                  "Unknown solver type \"" << solver << "\"." );
    }
    // ----------------------------------------------------------------------

    // Now the setting of inputs and outputs.
    EpetraExt::ModelEvaluator::InArgs inArgs = piro->createInArgs();
    RCP<Epetra_Vector> p1 =
      rcp(new Epetra_Vector(*(piro->get_p_init(0))));
    inArgs.set_p(0, p1);

    // Set output arguments to evalModel call.
    EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();

    // Now solve the problem and return the responses.
    const RCP<Teuchos::Time> piroSolveTime =
      Teuchos::TimeMonitor::getNewTimer("Piro total solve time");;
    {
    Teuchos::TimeMonitor tm(*piroSolveTime);
    piro->evalModel(inArgs, outArgs);
    }

    // Manually release LOCA stepper.
#ifdef HAVE_LOCA_ANASAZI
    if ( !glEigenSaver.is_null() )
      glEigenSaver->releaseLocaStepper();
#endif

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
