// @HEADER
//
//    Main continuation routine.
//    Copyright (C) 2012--2014  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
#include <string>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <NOX_Thyra_Vector.H>
#include <Piro_NOXSolver.hpp>
#include <Piro_LOCASolver.hpp>

#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Constant.hpp"
#include "nosh/ScalarField_ExplicitValues.hpp"
#include "nosh/ParameterMatrix_Keo.hpp"
#include "nosh/ParameterMatrix_DKeoDP.hpp"
#include "nosh/ParameterMatrix_Laplace.hpp"
#include "nosh/VectorField_ExplicitValues.hpp"
#include "nosh/VectorField_ConstantCurl.hpp"
#include "nosh/ModelEvaluator_Nls.hpp"
#include "nosh/Observer.hpp"
#include "nosh/SaveEigenData.hpp"
#include "nosh/CsvWriter.hpp"

#include "MyScalarField.hpp"

int main(int argc, char *argv[])
{
  // Create output stream. (Handy for multicore output.)
  auto out = Teuchos::VerboseObjectBase::getDefaultOStream();

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  auto comm = Teuchos::DefaultComm<int>::getComm();

  // Wrap the whole code in a big try-catch-statement.
  bool success = true;
  try {
    // =========================================================================
    // Handle command line arguments.
    // Boost::program_options is somewhat more complete here (e.g. you can
    // specify options without the "--" syntax), but it isn't less complicated
    // to use. Stick with Teuchos for now.
    Teuchos::CommandLineProcessor myClp;

    myClp.setDocString(
      "Numerical parameter continuation for nonlinear Schr\"odinger equations.\n"
    );

    std::string xmlInputPath = "";
    myClp.setOption("xml-input-file", &xmlInputPath,
                    "XML file containing the parameter list", true );

    // Print warning for unrecognized arguments and make sure to throw an
    // exception if something went wrong.
    //myClp.throwExceptions(false);
    //myClp.recogniseAllOptions ( true );

    // Finally, parse the command line.
    myClp.parse(argc, argv);

    // Retrieve Piro parameter list from given file.
    std::shared_ptr<Teuchos::ParameterList> piroParams(
        new Teuchos::ParameterList()
        );
    Teuchos::updateParametersFromXmlFile(
        xmlInputPath,
        Teuchos::rcp(piroParams).ptr()
        );
    // =======================================================================
    // Extract the location of input and output files.
    const Teuchos::ParameterList outputList =
      piroParams->sublist("Output", true);

    // Set default directory to be the directory of the XML file itself
    const std::string xmlDirectory =
      xmlInputPath.substr(0, xmlInputPath.find_last_of( "/" ) + 1);

    // By default, take the current directory.
    std::string prefix = "./";
    if (!xmlDirectory.empty())
      prefix = xmlDirectory + "/";

    const std::string outputDirectory = prefix;

    const std::string contFilePath =
      prefix + outputList.get<std::string>("Continuation data file name");

    Teuchos::ParameterList & inputDataList = piroParams->sublist("Input", true);

    const std::string inputExodusFile =
      prefix + inputDataList.get<std::string>("File");
    const int step = inputDataList.get<int>("Initial Psi Step");

    //const bool useBordering = piroParams->get<bool>("Bordering");
    // =======================================================================
    // Read the data from the file.
    auto mesh = std::make_shared<Nosh::StkMesh>(
        Teuchos::get_shared_ptr(comm),
        inputExodusFile,
        step
        );

    // Cast the data into something more accessible.
    auto psi = mesh->createComplexVector("psi");
    //psi->Random();

    // Set the output directory for later plotting with this.
    std::stringstream outputFile;
    outputFile << outputDirectory << "/solution.e";
    mesh->openOutputChannel(outputFile.str());

    // Create a parameter map from the initial parameter values.
    Teuchos::ParameterList initialParameterValues =
      piroParams->sublist("Initial parameter values", true);

    // Check if we need to interpret the time value stored in the file
    // as a parameter.
    const std::string & timeName =
      piroParams->get<std::string>("Interpret time as", "");
    if (!timeName.empty()) {
      initialParameterValues.set(timeName, mesh->getTime());
    }

    // Explicitly set the initial parameter value for this list.
    const std::string & paramName =
      piroParams->sublist( "LOCA" )
      .sublist( "Stepper" )
      .get<std::string>("Continuation Parameter");
    *out << "Setting the initial parameter value of \""
         << paramName << "\" to " << initialParameterValues.get<double>(paramName) << "." << std::endl;
    piroParams->sublist( "LOCA" )
    .sublist( "Stepper" )
    .set("Initial Value", initialParameterValues.get<double>(paramName));

    // Set the thickness field.
    auto thickness = std::make_shared<Nosh::ScalarField::Constant>(*mesh, 1.0);

    // Some alternatives for the positive-definite operator.
    // (a) -\Delta (Laplace operator with Neumann boundary)
    //const std::shared_ptr<Nosh::ParameterMatrix::Virtual> matrixBuilder =
    //  rcp(new Nosh::ParameterMatrix::Laplace(mesh, thickness));

    // (b) (-i\nabla-A)^2 (Kinetic energy of a particle in magnetic field)
    // (b1) 'A' explicitly given in file.
    const double mu = initialParameterValues.get<double>("mu");
    auto mvp = std::make_shared<Nosh::VectorField::ExplicitValues>(*mesh, "A", mu);

    //const std::shared_ptr<Nosh::ParameterMatrix::Virtual> keoBuilder(
    //    new Nosh::ParameterMatrix::Keo(mesh, thickness, mvp)
    //    );
    //const std::shared_ptr<Nosh::ParameterMatrix::Virtual> DKeoDPBuilder(
    //    new Nosh::ParameterMatrix::DKeoDP(mesh, thickness, mvp, "mu")
    //    );

    // (b2) 'A' analytically given (here with constant curl).
    //      Optionally add a rotation axis u. This is important
    //      if continuation happens as a rotation of the vector
    //      field around an axis.
    //const std::shared_ptr<DoubleVector> b = rcp(new DoubleVector(3));
    //std::shared_ptr<Teuchos::SerialDenseVector<int,double> > u = Teuchos::null;
    //if ( piroParams->isSublist("Rotation vector") )
    //{
    //    u = rcp(new Teuchos::SerialDenseVector<int,double>(3));
    //    Teuchos::ParameterList & rotationVectorList =
    //        piroParams->sublist( "Rotation vector", false );
    //    (*u)[0] = rotationVectorList.get<double>("x");
    //    (*u)[1] = rotationVectorList.get<double>("y");
    //    (*u)[2] = rotationVectorList.get<double>("z");
    //}
    //std::shared_ptr<Nosh::VectorField::Virtual> mvp =
    //  rcp(new Nosh::VectorField::ConstantCurl(mesh, b, u));
    //const std::shared_ptr<Nosh::ParameterMatrix::Virtual> matrixBuilder =
    //  rcp(new Nosh::ParameterMatrix::Keo(mesh, thickness, mvp));
    // (b3) 'A' analytically given in a class you write yourself, derived
    //      from Nosh::ParameterMatrix::Virtual.
    // [...]
    //
    // Setup the scalar potential V.
    // (a) A constant potential.
    //std::shared_ptr<Nosh::ScalarField::Virtual> sp =
    //rcp(new Nosh::ScalarField::Constant(*mesh, -1.0));
    //const double T = initialParameterValues.get<double>("T");
    // (b) With explicit values.
    //std::shared_ptr<Nosh::ScalarField::Virtual> sp =
    //rcp(new Nosh::ScalarField::ExplicitValues(*mesh, "V"));
    // (c) One you built yourself by deriving from Nosh::ScalarField::Virtual.
    auto sp = std::make_shared<MyScalarField>(mesh);


    const double g = initialParameterValues.get<double>("g");
    // Finally, create the model evaluator.
    // This is the most important object in the whole stack.
    auto modelEvaluator = std::make_shared<Nosh::ModelEvaluator::Nls>(
        mesh,
        mvp,
        sp,
        g,
        thickness,
        psi,
        "mu"
        );

    // Build the Piro model evaluator. It's used to hook up with
    // several different backends (NOX, LOCA, Rhythmos,...).
    std::shared_ptr<Thyra::ModelEvaluator<double>> piro;

    // Declare the eigensaver; it will be used only for LOCA solvers, though.
    std::shared_ptr<Nosh::SaveEigenData> glEigenSaver;

    // Switch by solver type.
    std::string & solver = piroParams->get<std::string>("Piro Solver");

    if (solver == "NOX") {
      auto observer = std::make_shared<Nosh::Observer>(modelEvaluator);

      piro = std::make_shared<Piro::NOXSolver<double>>(
            Teuchos::rcp(piroParams),
            Teuchos::rcp(modelEvaluator),
            Teuchos::rcp(observer)
            );
    } else if (solver == "LOCA") {
      auto observer = std::make_shared<Nosh::Observer>(
          modelEvaluator,
          contFilePath,
          piroParams->sublist("LOCA")
          .sublist("Stepper")
          .get<std::string>("Continuation Parameter")
          );

      // Setup eigen saver.
#ifdef HAVE_LOCA_ANASAZI
      bool computeEigenvalues = piroParams->sublist( "LOCA" )
                                .sublist( "Stepper" )
                                .get<bool>("Compute Eigenvalues");
      if (computeEigenvalues) {
        Teuchos::ParameterList & eigenList = piroParams->sublist("LOCA")
                                             .sublist("Stepper")
                                             .sublist("Eigensolver");
        std::string eigenvaluesFilePath =
          xmlDirectory + "/" + outputList.get<std::string> ( "Eigenvalues file name" );

        glEigenSaver = std::make_shared<Nosh::SaveEigenData>(
            eigenList,
            modelEvaluator,
            eigenvaluesFilePath
            );

        std::shared_ptr<LOCA::SaveEigenData::AbstractStrategy>
          glSaveEigenDataStrategy = glEigenSaver;
        eigenList.set("Save Eigen Data Method",
                      "User-Defined");
        eigenList.set("User-Defined Save Eigen Data Name",
                      "glSaveEigenDataStrategy");
        eigenList.set("glSaveEigenDataStrategy",
                      glSaveEigenDataStrategy);
      }
#endif
      // Get the solver.
      std::shared_ptr<Piro::LOCASolver<double>> piroLOCASolver(
          new Piro::LOCASolver<double>(
            Teuchos::rcp(piroParams),
            Teuchos::rcp(modelEvaluator),
            Teuchos::null
            //Teuchos::rcp(observer)
            )
          );

//      // Get stepper and inject it into the eigensaver.
//      std::shared_ptr<LOCA::Stepper> stepper = Teuchos::get_shared_ptr(
//          piroLOCASolver->getLOCAStepperNonConst()
//          );
//#ifdef HAVE_LOCA_ANASAZI
//      if (computeEigenvalues)
//        glEigenSaver->setLocaStepper(stepper);
//#endif
      piro = piroLOCASolver;
    }
#if 0
    else if ( solver == "Turning Point" ) {
      std::shared_ptr<Nosh::Observer> observer;

      Teuchos::ParameterList & bifList =
        piroParams->sublist("LOCA").sublist("Bifurcation");

      // Fetch the (approximate) null state.
      auto nullstateZ = mesh->createVector("null");

      // Set the length normalization vector to be the initial null vector.
      TEUCHOS_ASSERT(nullstateZ);
      auto lengthNormVec = Teuchos::rcp(new NOX::Thyra::Vector(*nullstateZ));
      //lengthNormVec->init(1.0);
      bifList.set("Length Normalization Vector", lengthNormVec);

      // Set the initial null vector.
      auto initialNullAbstractVec =
        Teuchos::rcp(new NOX::Thyra::Vector(*nullstateZ));
      // initialNullAbstractVec->init(1.0);
      bifList.set("Initial Null Vector", initialNullAbstractVec);

      piro = std::make_shared<Piro::LOCASolver<double>>(
            Teuchos::rcp(piroParams),
            Teuchos::rcp(modelEvaluator),
            Teuchos::null
            //Teuchos::rcp(observer)
            );
    }
#endif
    else {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          "Unknown solver type \"" << solver << "\"."
          );
    }
    // ----------------------------------------------------------------------

    // Now the setting of inputs and outputs.
    Thyra::ModelEvaluatorBase::InArgs<double> inArgs = piro->createInArgs();
    inArgs.set_p(
        0,
        piro->getNominalValues().get_p(0)
        );

    // Set output arguments to evalModel call.
    Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = piro->createOutArgs();

    // Now solve the problem and return the responses.
    const Teuchos::RCP<Teuchos::Time> piroSolveTime =
      Teuchos::TimeMonitor::getNewTimer("Piro total solve time");;
    {
      Teuchos::TimeMonitor tm(*piroSolveTime);
      piro->evalModel(inArgs, outArgs);
    }

    // Manually release LOCA stepper.
#ifdef HAVE_LOCA_ANASAZI
    if (glEigenSaver)
      glEigenSaver->releaseLocaStepper();
#endif

    // Print timing data.
    Teuchos::TimeMonitor::summarize();
  } catch (Teuchos::CommandLineProcessor::HelpPrinted) {
  } catch (Teuchos::CommandLineProcessor::ParseError) {
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
