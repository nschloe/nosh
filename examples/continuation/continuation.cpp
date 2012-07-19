// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Vector.h>

#include <Piro_Epetra_NOXSolver.hpp>
#include <Piro_Epetra_LOCASolver.hpp>

#include "Nosh_StkMesh.hpp"
#include "Nosh_StkMeshReader.hpp"
#include "Nosh_State.hpp"
#include "Nosh_ScalarPotential_Constant.hpp"
#include "Nosh_MatrixBuilder_Keo.hpp"
#include "Nosh_MagneticVectorPotential_ExplicitValues.hpp"
#include "Nosh_MagneticVectorPotential_ConstantField.hpp"
#include "Nosh_ModelEvaluator.hpp"
#include "Nosh_NoxObserver.hpp"
#include "Nosh_SaveEigenData.hpp"
#include "Nosh_CsvWriter.hpp"

#include <Teuchos_TimeMonitor.hpp>

// =============================================================================
std::string
extractDirectory(const std::string& path)
{
  // Extract the directory from a string, e.g.,
  // "/usr/bin/" from "/usr/bin/gcc".
  return path.substr(0, path.find_last_of( "/" ) + 1);
}
// =============================================================================
int main(int argc, char *argv[])
{
  // Initialize MPI.
#ifdef HAVE_MPI
  MPI_Init( &argc, &argv );
#endif

  // Create a communicator for Epetra objects.
#ifdef HAVE_MPI
  Teuchos::RCP<Epetra_MpiComm> eComm =
    Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
  Teuchos::RCP<Epetra_SerialComm> eComm =
    Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

  // Create output stream. (Handy for multicore output.)
  const Teuchos::RCP<Teuchos::FancyOStream> out =
      Teuchos::VerboseObjectBase::getDefaultOStream();

  bool success = true;
  try
  {
    // ===========================================================================
    // Handle command line arguments.
    Teuchos::CommandLineProcessor My_CLP;

    My_CLP.setDocString(
        "This program solves nonlinear Schrödinger equations with a Piro interface.\n"
    );

    std::string xmlInputFileName = "";
    My_CLP.setOption ("xml-input-file", &xmlInputFileName,
                      "XML file containing the parameter list", true );

    // Print warning for unrecognized arguments and make sure to throw an
    // exception if something went wrong.
    My_CLP.recogniseAllOptions ( true );
    My_CLP.throwExceptions(true);

    // Finally, parse the command line.
    My_CLP.parse(argc, argv);

    // Retrieve Piro parmeter list from given file.
    Teuchos::RCP<Teuchos::ParameterList> piroParams =
        Teuchos::rcp(new Teuchos::ParameterList());
    Teuchos::updateParametersFromXmlFile(xmlInputFileName,
                                         piroParams.ptr());
    // =======================================================================
    // Extract the location of input and output files.
    const Teuchos::ParameterList outputList = piroParams->sublist("Output", true);

    // Set default directory to be the directory of the XML file itself
    const std::string xmlDirectory = extractDirectory( xmlInputFileName );

    const std::string & outputDirectory = xmlDirectory;

    const std::string contFilePath = xmlDirectory + "/"
                                   + outputList.get<std::string>( "Continuation data file name" );

    Teuchos::ParameterList & inputDataList =
      piroParams->sublist ( "Input", true );
    const std::string constantDataFile = xmlDirectory + "/"
                                       + inputDataList.get<std::string>( "Constant data" );
    const std::string initialPsiFile = xmlDirectory + "/"
                                     + inputDataList.sublist("Initial psi", true)
                                       .get<std::string>( "File" );
    const int step = inputDataList.sublist("Initial psi", true).get<int>( "Step" );
    // =======================================================================
    // Get the initial parameter values.
    Teuchos::ParameterList initialParameterValues =
      piroParams->sublist("Initial parameter values", true);

    // Read the the non-transient data from the first time step.
    Teuchos::ParameterList data;
    Nosh::StkMeshRead(*eComm, constantDataFile, 0, data);
    // Cast the data into something more accessible.
    const Teuchos::RCP<Nosh::StkMesh> & mesh =
      data.get<Teuchos::RCP<Nosh::StkMesh> >( "mesh" );
    const Teuchos::RCP<const Epetra_MultiVector> & mvpValues =
      data.get("A", Teuchos::RCP<const Epetra_MultiVector>() );
    const Teuchos::RCP<Epetra_Vector> & potentialValues =
      data.get("V", Teuchos::RCP<Epetra_Vector>());
    const Teuchos::RCP<Epetra_Vector> & thickness =
      data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );

    // Read PSI from a given time step.
    Teuchos::RCP<Epetra_Vector> psi;
    if (constantDataFile.compare(initialPsiFile) && step == 0)
      psi = data.get("psi", Teuchos::RCP<Epetra_Vector>() );
    else
    {
      Teuchos::ParameterList tData;
      Nosh::StkMeshRead(*eComm, initialPsiFile, step, tData);
      psi = tData.get("psi", Teuchos::RCP<Epetra_Vector>() );
    }

    // Set the output directory for later plotting with this.
    mesh->openOutputChannel(outputDirectory, "solution");

    // Create the initial state from psi.
    TEUCHOS_ASSERT( !psi.is_null() );
    Teuchos::RCP<Nosh::State> state =
      Teuchos::rcp(new Nosh::State(*psi, mesh));

    // Setup the energy operator, here: (-i\nabla-A)^2.
    // Choose between several given MVPs or build your own by
    // deriving from Nosh::MagneticVectorPotential::Virtual.
    const double initMu = initialParameterValues.get<double>("mu", 0.0);
    Teuchos::RCP<Nosh::MagneticVectorPotential::Virtual> mvp =
      Teuchos::rcp(new Nosh::MagneticVectorPotential::ExplicitValues(mesh, mvpValues, initMu));
    //const Teuchos::RCP<DoubleVector> b = Teuchos::rcp(new DoubleVector(3));
    //Teuchos::RCP<Nosh::MagneticVectorPotential::Virtual> mvp =
    //  Teuchos::rcp(new Nosh::MagneticVectorPotential::ConstantField(mesh, b));
    const Teuchos::RCP<Nosh::MatrixBuilder::Virtual> matrixBuilder =
        Teuchos::rcp(new Nosh::MatrixBuilder::Keo(mesh, thickness, mvp));

    // Alternative: Analytically given MVP. This one can also be rotated in space.
    //const double theta = initialParameterValues.get<double>("theta", 0.0);
    // Get the rotation vector.
    // This is important if continuation happens as a rotation of the
    // vector field around an axis.
    //Teuchos::RCP<Teuchos::SerialDenseVector<int,double> > u = Teuchos::null;
    //if ( piroParams->isSublist("Rotation vector") )
    //{
    //    u = Teuchos::rcp(new Teuchos::SerialDenseVector<int,double>(3) );
    //    Teuchos::ParameterList & rotationVectorList =
    //        piroParams->sublist( "Rotation vector", false );
    //    (*u)[0] = rotationVectorList.get<double>("x");
    //    (*u)[1] = rotationVectorList.get<double>("y");
    //    (*u)[2] = rotationVectorList.get<double>("z");
    //}
    //Teuchos::RCP<Nosh::MagneticVectorPotential::Virtual> mvp =
    //  Teuchos::rcp ( new Nosh::MagneticVectorPotential::ConstantField(mesh, mvpValues, mu, theta, u));

    // Setup the scalar potential V.
    // Use this or build your own by deriving from Nosh::ScalarPotential::Virtual.
    const double T = initialParameterValues.get<double>("T", 0.0);
    Teuchos::RCP<Nosh::ScalarPotential::Virtual> sp =
      Teuchos::rcp(new Nosh::ScalarPotential::Constant(-1.0));

    // Finally, create the model evaluator.
    // This is the most important object in the whole stack.
    const double g = initialParameterValues.get<double>("g");
    Teuchos::RCP<Nosh::ModelEvaluator> nlsModel =
      Teuchos::rcp(new Nosh::ModelEvaluator(mesh, matrixBuilder, sp, g, thickness, psi));

    // Build the Piro model evaluator. It's used to hook up with
    // several different backends (NOX, LOCA, Rhythmos,...).
    Teuchos::RCP<EpetraExt::ModelEvaluator> piro;

    // Declare the eigensaver; it will be used only for LOCA solvers, though.
    Teuchos::RCP<Nosh::SaveEigenData> glEigenSaver;

    // Switch by solver type.
    std::string & solver = piroParams->get("Piro Solver", "");
    // ----------------------------------------------------------------------
    if (solver == "NOX")
    {
      Teuchos::RCP<Nosh::NoxObserver> observer =
        Teuchos::rcp(new Nosh::NoxObserver(nlsModel,
                                            contFilePath,
                                            Nosh::NoxObserver::OBSERVER_TYPE_NEWTON
                                            ));

      piro = Teuchos::rcp(new Piro::Epetra::NOXSolver(piroParams,
                                                      nlsModel,
                                                      observer));
    }
    // ----------------------------------------------------------------------
    else if (solver == "LOCA")
    {
      Teuchos::RCP<Nosh::NoxObserver> observer =
        Teuchos::rcp(new Nosh::NoxObserver(nlsModel,
                                            contFilePath,
                                            Nosh::NoxObserver::OBSERVER_TYPE_CONTINUATION
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

        Teuchos::RCP<Nosh::CsvWriter> eigenCsvWriter =
          Teuchos::rcp( new Nosh::CsvWriter( eigenvaluesFilePath ) );

        glEigenSaver =
          Teuchos::RCP<Nosh::SaveEigenData>(new Nosh::SaveEigenData(eigenList,
                                                                      nlsModel,
                                                                      eigenCsvWriter));

        Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy =
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
      Teuchos::RCP<Piro::Epetra::LOCASolver> piroLOCASolver =
        Teuchos::rcp( new Piro::Epetra::LOCASolver(piroParams, nlsModel, observer));

      // Get stepper and inject it into the eigensaver.
      Teuchos::RCP<LOCA::Stepper> stepper = piroLOCASolver->getLOCAStepperNonConst();
 #ifdef HAVE_LOCA_ANASAZI
      if (computeEigenvalues)
        glEigenSaver->setLocaStepper( stepper );
 #endif
      piro = piroLOCASolver;
    }
    // ----------------------------------------------------------------------
    else if ( solver == "Turning Point" )
    {
      Teuchos::RCP<Nosh::NoxObserver> observer = Teuchos::null;

      // Read the data from the file.
      std::string nullstateFilePath = xmlDirectory + "/" + inputDataList.get<std::string> ( "Null state" );
      Teuchos::ParameterList nullstateData;
      Nosh::StkMeshRead( *eComm, nullstateFilePath, 0, nullstateData );

      // Cast the data into something more accessible.
      Teuchos::RCP<Nosh::StkMesh> & nullstateMesh = nullstateData.get( "mesh", Teuchos::RCP<Nosh::StkMesh>() );
      Teuchos::RCP<Epetra_Vector>  & nullstateZ = nullstateData.get( "psi", Teuchos::RCP<Epetra_Vector>() );

      TEUCHOS_ASSERT( !nullstateZ.is_null() );
      Teuchos::RCP<Nosh::State> nullstate =
          Teuchos::rcp( new Nosh::State( *nullstateZ, nullstateMesh ) );

      Teuchos::ParameterList & bifList =
          piroParams->sublist ( "LOCA" ).sublist ( "Bifurcation" );

      // Set the length normalization vector to be the initial null vector.
      Teuchos::RCP<NOX::Abstract::Vector> lengthNormVec =
          Teuchos::rcp ( new NOX::Epetra::Vector ( *(nullstate->getPsi()) ) );
      //lengthNormVec->init(1.0);
      bifList.set ( "Length Normalization Vector", lengthNormVec );

      // Set the initial null vector.
      Teuchos::RCP<NOX::Abstract::Vector> initialNullAbstractVec =
          Teuchos::rcp(new NOX::Epetra::Vector(*(nullstate->getPsi())));
      // initialNullAbstractVec->init(1.0);
      bifList.set ( "Initial Null Vector", initialNullAbstractVec );

      piro = Teuchos::rcp(new Piro::Epetra::LOCASolver( piroParams,
                                                        nlsModel,
                                                        observer ));
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
    Teuchos::RCP<Epetra_Vector> p1 =
        Teuchos::rcp(new Epetra_Vector(*(piro->get_p_init(0))));
    inArgs.set_p(0, p1);

    // Set output arguments to evalModel call
    EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();

    // Now solve the problem and return the responses.
    const Teuchos::RCP<Teuchos::Time> piroSolveTime =
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
