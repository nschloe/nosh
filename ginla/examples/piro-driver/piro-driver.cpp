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

#include <LOCA_StatusTest_Combo.H>

#include <Piro_Epetra_NOXSolver.hpp>
#include <Piro_Epetra_LOCASolver.hpp>

#include "Ginla_StkMesh.hpp"
#include "Ginla_StkMeshReader.hpp"
#include "Ginla_State.hpp"
#include "Ginla_MagneticVectorPotential.hpp"
#include "Ginla_ModelEvaluator.hpp"
#include "Ginla_NoxObserver.hpp"
#include "Ginla_SaveEigenData.hpp"
#include "Ginla_StateWriter.hpp"

#include <Teuchos_TimeMonitor.hpp>

// =============================================================================
std::string
extractDirectory( const std::string& path )
{
    return path.substr( 0, path.find_last_of( "/" ) +1 );
}
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
          "This program solves the Ginzburg--Landau problem with a Piro interface.\n"
      );

      std::string xmlInputFileName = "";
      My_CLP.setOption ( "xml-input-file", &xmlInputFileName,
                         "XML file containing the parameter list", true );

      // print warning for unrecognized arguments
      My_CLP.recogniseAllOptions ( true );

      // do throw exceptions
      My_CLP.throwExceptions ( true );

      // finally, parse the command line
      Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
      My_CLP.parse ( argc, argv );

      Teuchos::RCP<Teuchos::ParameterList> piroParams =
          Teuchos::rcp ( new Teuchos::ParameterList );
      if ( eComm->MyPID() == 0 )
          std::cout << "Reading parameter list from \"" << xmlInputFileName << "\"."
                    << std::endl;

      Teuchos::updateParametersFromXmlFile( xmlInputFileName, piroParams.get() );

      // =======================================================================
      // extract data of the parameter list
      Teuchos::ParameterList outputList = piroParams->sublist ( "Output", true );

      // set default directory to be the directory of the XML file itself
      std::string xmlDirectory = extractDirectory( xmlInputFileName );

      std::string & outputDirectory = xmlDirectory;

      std::string contFilePath = xmlDirectory + "/"
                               + outputList.get<std::string> ( "Continuation data file name" );

      Teuchos::ParameterList initialGuessList;
      initialGuessList = piroParams->sublist ( "Initial guess", true );
      std::string inputFilePath = xmlDirectory + "/"
                                + initialGuessList.get<std::string> ( "State" );
      // =======================================================================
      // Read the data from the file.
      Teuchos::ParameterList data;
      Ginla::StkMeshRead( *eComm, inputFilePath, data );

      // Cast the data into something more accessible.
      Teuchos::RCP<Ginla::StkMesh>     & mesh = data.get<Teuchos::RCP<Ginla::StkMesh> >( "mesh" );
      Teuchos::RCP<Epetra_Vector>      & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
      Teuchos::RCP<const Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<const Epetra_MultiVector>() );
      Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );

      // set the output directory for later plotting with this
      mesh->setOutputFile( outputDirectory, "solution" );

      // create the state
      TEUCHOS_ASSERT( !z.is_null() );
      Teuchos::RCP<Ginla::State> state =
          Teuchos::rcp( new Ginla::State( *z, mesh ) );

      // possibly override the parameters
      Teuchos::ParameterList problemParameters = Teuchos::ParameterList();
      Teuchos::ParameterList & OverrideParamsList =
          piroParams->sublist( "Override parameter list", true );
      problemParameters.setParameters( OverrideParamsList );

      // Get the rotation vector.
      Teuchos::RCP<Teuchos::SerialDenseVector<int,double> > u = Teuchos::null;
      if ( piroParams->isSublist("Rotation vector") )
      {
          u = Teuchos::rcp(new Teuchos::SerialDenseVector<int,double>(3) );
          Teuchos::ParameterList & rotationVectorList =
              piroParams->sublist( "Rotation vector", false );
          (*u)[0] = rotationVectorList.get<double>("x");
          (*u)[1] = rotationVectorList.get<double>("y");
          (*u)[2] = rotationVectorList.get<double>("z");
      }

      double mu = problemParameters.get<double>( "mu", 0.0 );
      double theta = problemParameters.get<double>( "theta", 0.0 );
      Teuchos::RCP<Ginla::MagneticVectorPotential> mvp =
              Teuchos::rcp ( new Ginla::MagneticVectorPotential ( mesh, mvpValues, mu, theta, u ) );

      // create the mode evaluator
      Teuchos::RCP<Ginla::ModelEvaluator> glModel =
              Teuchos::rcp( new Ginla::ModelEvaluator( mesh,
                                                       problemParameters,
                                                       thickness,
                                                       mvp,
                                                       z
                                                     )
                          );

      Teuchos::RCP<Ginla::NoxObserver> observer;

      Teuchos::RCP<Ginla::StatsWriter> statsWriter =
          Teuchos::rcp( new Ginla::StatsWriter( contFilePath ) );

      // warn if initial value was given twice
      std::string contParam = piroParams->sublist ( "LOCA" )
          .sublist ( "Stepper" ).get<std::string> ( "Continuation Parameter" );
      if ( problemParameters.isParameter ( contParam ) )
          *out << "Warning: Continuation parameter \""
               << contParam
               << "\" explicitly given. Initial value will be overwritten by "
               << "'LOCA->Stepper->Initial Value', though."
               << std::endl;

      // Use these two objects to construct a Piro solved application
      //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
      Teuchos::RCP<EpetraExt::ModelEvaluator> piro;

      // Declare the eigensaver; it will be used only for LOCA solvers, though.
      Teuchos::RCP<Ginla::SaveEigenData> glEigenSaver;

      // Setup the data output.
      int maxLocaSteps = piroParams->sublist ( "LOCA" )
                                    .sublist ( "Stepper" )
                                    .get<int> ( "Max Steps" );
      Teuchos::RCP<Ginla::StateWriter> stateWriter =
          Teuchos::rcp( new Ginla::StateWriter( outputDirectory,
                                                "solution"
                                              )
                      );

      // switch by solver type
      std::string & solver = piroParams->get( "Piro Solver", "" );
      // ----------------------------------------------------------------------
      if ( solver == "NOX" )
      {
          observer = Teuchos::rcp( new Ginla::NoxObserver( glModel,
                                                           Ginla::NoxObserver::OBSERVER_TYPE_NEWTON
                                                         )
                                 );
          observer->setStatisticsWriter( statsWriter );

          piro = Teuchos::rcp(new Piro::Epetra::NOXSolver( piroParams,
                                                           glModel,
                                                           observer ));
      }
      // ----------------------------------------------------------------------
      else if ( solver == "LOCA" )
      {
          observer = Teuchos::rcp( new Ginla::NoxObserver( glModel,
                                                           Ginla::NoxObserver::OBSERVER_TYPE_CONTINUATION
                                                         )
                                 );
          observer->setStatisticsWriter( statsWriter );

//          Teuchos::RCP<LOCA::StatusTest::Combo> locaTest =
//              Teuchos::rcp( new LOCA::StatusTest::Combo( LOCA::StatusTest::Combo::OR ) );

          // setup eingen saver
#ifdef HAVE_LOCA_ANASAZI
           Teuchos::ParameterList & eigenList = piroParams->sublist( "LOCA" )
                                                           .sublist( "Stepper" )
                                                           .sublist( "Eigensolver" );
           std::string eigenvaluesFilePath = xmlDirectory
                                           + "/"
                                           + outputList.get<std::string> ( "Eigenvalues file name" );

           Teuchos::RCP<Ginla::StatsWriter> eigenStatsWriter =
               Teuchos::rcp( new Ginla::StatsWriter( eigenvaluesFilePath ) );

           // initialize the stability change test with a pointer to the eigenvalue information
           int stabilityChangeTreshold = 1; // stop when the stability changes by multiplicity 1
           Teuchos::RCP<const Teuchos::ParameterList> eigendataList = eigenStatsWriter->getList();

           //Teuchos::RCP<LOCA::StatusTest::Abstract> stabilityChangeTest =
           //    Teuchos::rcp( new Ginla::StatusTest::StabilityChange( eigendataList,
           //                                                          stabilityChangeTreshold ) );
           //locaTest->addStatusTest( stabilityChangeTest );

           glEigenSaver = Teuchos::RCP<Ginla::SaveEigenData> ( new Ginla::SaveEigenData ( eigenList,
                                                                                          glModel,
                                                                                          eigenStatsWriter ) );

           Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy = glEigenSaver;
           eigenList.set ( "Save Eigen Data Method", "User-Defined" );
           eigenList.set ( "User-Defined Save Eigen Data Name", "glSaveEigenDataStrategy" );
           eigenList.set ( "glSaveEigenDataStrategy", glSaveEigenDataStrategy );
#endif

//           Teuchos::RCP<LOCA::StatusTest::Abstract> freeEnergyTest =
//               Teuchos::rcp( new Ginla::StatusTest::Energy( glModel, 0.0 ) );
//           locaTest->addStatusTest( freeEnergyTest );
//
//           // feed in restart predictor vector
//           Teuchos::ParameterList & predictorList = piroParams->sublist ( "LOCA" )
//                                                               .sublist ( "Predictor" );
//           if ( predictorList.get<std::string>("Method") == "Secant" )
//           {
//               Teuchos::ParameterList & fspList = predictorList.sublist ( "First Step Predictor" );
//               if ( fspList.get<std::string>("Method") == "Restart" )
//               {
//                   // Get restart vector.
//                   // read the predictor state
//                   Teuchos::ParameterList           voidParameters;
//                   Teuchos::RCP<ComplexMultiVector> z = Teuchos::null;
//                   Teuchos::RCP<VMesh::Mesh>    mesh = Teuchos::null;
//                   VMesh::read( Comm,
//                                    getAbsolutePath(initialGuessList.get<std::string> ( "Predictor" ), xmlPath),
//                                    z,
//                                    mesh,
//                                    voidParameters );
//                   Teuchos::RCP<Ginla::FVM::State>  predictorState =
//                       Teuchos::rcp( new Ginla::FVM::State( z, mesh ) );
//
//                   // transform to system vector
//                   Teuchos::RCP<Epetra_Vector> predictorV = glModel->createSystemVector( *predictorState );
//
//                   // Create the LOCA::MultiContinuation::ExtendedVector;
//                   // that's the predictor vector with a prediction for the parameter appended
//                   // -- 0 in the case of a pitchfork.
//                   NOX::Epetra::Vector vx( predictorV );
//
//                   // Note that the globalData element is here the null pointer. It would be work-aroundish to
//                   // artificially create one here, and after all, all that's it's used for is error and warning
//                   // message printing.
//                   Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> restartVector =
//                       Teuchos::rcp( new LOCA::MultiContinuation::ExtendedVector( Teuchos::null,
//                                                                                  vx,
//                                                                                  1 )
//                                   );
//                   double vp = 0.0;
//                   restartVector->setScalar(0, vp);
//
//                   fspList.set( "Restart Vector", restartVector );
//
// //                   // Make only one LOCA step when branch switching.
// //                   bool return_failed_on_max_steps = false;
// //                   extraTest = Teuchos::rcp( new Ginla::StatusTest::MaxAcceptedSteps( 2,
// //                                                                                      return_failed_on_max_steps
// //                                                                                     ) );
//               }
//           }
//
//           double lowerLimit = stepperList.get<double> ( "Min Value" );
//           double upperLimit = stepperList.get<double> ( "Max Value" );
//           Teuchos::RCP<Ginla::StatusTest::ParameterLimits> paramLimitsTest =
//               Teuchos::rcp( new Ginla::StatusTest::ParameterLimits( lowerLimit, upperLimit, false ) );
//
//           locaTest->addStatusTest( paramLimitsTest );

          // fetch the stepper
          Teuchos::RCP<Piro::Epetra::LOCASolver> piroLOCASolver =
              Teuchos::rcp( new Piro::Epetra::LOCASolver( piroParams,
                                                          glModel,
                                                          observer//,
//                                                           Teuchos::null,
//                                                           locaTest
                                                        ) );

           // get stepper and inject it into the eigensaver
           Teuchos::RCP<LOCA::Stepper> stepper = piroLOCASolver->getLOCAStepperNonConst();
 #ifdef HAVE_LOCA_ANASAZI
           glEigenSaver->setLocaStepper ( stepper );
 #endif
          piro = piroLOCASolver;
      }
      // ----------------------------------------------------------------------
      else if ( solver == "Turning Point" )
      {
          // TODO make sure the turning point continuation doesn't technically fail by default

          observer = Teuchos::null;
//          observer = Teuchos::rcp( new Ginla::NoxObserver( stateWriter,
//                                                               glModel,
//                                                               Ginla::NoxObserver::OBSERVER_TYPE_TURNING_POINT ) );
//           observer->setStatisticsWriter( statsWriter, glOperator );

          // get the initial null state file
          initialGuessList = piroParams->sublist ( "Initial guess", true );

          // Read the data from the file.
          std::string nullstateFilePath = xmlDirectory + "/" + initialGuessList.get<std::string> ( "Null state" );
          Teuchos::ParameterList data;
          Ginla::StkMeshRead( *eComm, nullstateFilePath, data );

          // Cast the data into something more accessible.
          Teuchos::RCP<Ginla::StkMesh> & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );
          Teuchos::RCP<Epetra_Vector>  & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );

          TEUCHOS_ASSERT( !z.is_null() );
          Teuchos::RCP<Ginla::State> nullstate =
              Teuchos::rcp( new Ginla::State( *z, mesh ) );

          Teuchos::ParameterList & bifList =
              piroParams->sublist ( "LOCA" ).sublist ( "Bifurcation" );

          // set the length normalization vector to be the initial null vector
          Teuchos::RCP<NOX::Abstract::Vector> lengthNormVec =
              Teuchos::rcp ( new NOX::Epetra::Vector ( *(nullstate->getPsi()) ) );
  //         lengthNormVec->init(1.0);
          bifList.set ( "Length Normalization Vector", lengthNormVec );

          // set the initial null vector
          Teuchos::RCP<NOX::Abstract::Vector> initialNullAbstractVec =
              Teuchos::rcp ( new NOX::Epetra::Vector ( *(nullstate->getPsi()) ) );
      //     initialNullAbstractVec->init(1.0);
          bifList.set ( "Initial Null Vector", initialNullAbstractVec );

          piro = Teuchos::rcp(new Piro::Epetra::LOCASolver( piroParams,
                                                            glModel,
                                                            observer ));
      }
      // ----------------------------------------------------------------------
      else
      {
          TEST_FOR_EXCEPT_MSG( true,
                               "Unknown solver type \"" << solver << "\"." );
      }
      // ----------------------------------------------------------------------

      // Now the (somewhat cumbersome) setting of inputs and outputs
      EpetraExt::ModelEvaluator::InArgs inArgs = piro->createInArgs();
      int num_p = inArgs.Np();     // Number of *vectors* of parameters
      Teuchos::RCP<Epetra_Vector> p1 =
          Teuchos::rcp(new Epetra_Vector(*(piro->get_p_init(0))));
      inArgs.set_p( 0, p1 );

      // Set output arguments to evalModel call
      EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();

      // Now, solve the problem and return the responses
      const Teuchos::RCP<Teuchos::Time> piroSolveTime =
          Teuchos::TimeMonitor::getNewTimer("Piro total solve time");;
      {
      Teuchos::TimeMonitor tm(*piroSolveTime);
      piro->evalModel(inArgs, outArgs);
      }
      // Make sure it finsihsed without error.
      //TEUCHOS_ASSERT( !outArgs.isFailed() );

      // manually release LOCA stepper
#ifdef HAVE_LOCA_ANASAZI
      if ( !glEigenSaver.is_null() )
          glEigenSaver->releaseLocaStepper();
#endif

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
