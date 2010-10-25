#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Piro_Epetra_NOXSolver.hpp>
#include <Piro_Epetra_LOCASolver.hpp>

#include <boost/filesystem.hpp>

#include "Ginla_EpetraFVM_ModelEvaluator.h"
#include "Ginla_EpetraFVM_State.h"
#include "VIO_EpetraMesh_Reader.h"

// #include "Ginla_IO_SaveNewtonData.h"
#include "Ginla_IO_SaveEigenData.h"
#include "Ginla_IO_NoxObserver.h"

#include "Ginla_MagneticVectorPotential_X.h"
#include "Ginla_MagneticVectorPotential_Y.h"
#include "Ginla_MagneticVectorPotential_Z.h"
// #include "Ginla_MagneticVectorPotential_ZSquareSymmetric.h"

#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_StatsWriter.h"

#include "Ginla_StatusTest_MaxAcceptedSteps.h"
#include "Ginla_StatusTest_Energy.h"
#include "Ginla_StatusTest_ParameterLimits.h"
#include "Ginla_StatusTest_StabilityChange.h"
#include "LOCA_StatusTest_Combo.H"

// =============================================================================
// declarations (definitions below)
std::string
getAbsolutePath(       boost::filesystem::path   filepath,
                 const boost::filesystem::path & xmlPath
               )
{
    if ( !filepath.empty() && filepath.root_directory().empty() ) // filepath is a relative path
        filepath = xmlPath / filepath;
    TEUCHOS_ASSERT( !filepath.empty() );

    return filepath.string();
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

    int status = 0;

    try
    {
      // ===========================================================================
      // handle command line arguments
      Teuchos::CommandLineProcessor My_CLP;

      My_CLP.setDocString (
          "This program solves the Ginzburg--Landau problem with a Piro LOCA interface.\n"
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

      Teuchos::updateParametersFromXmlFile ( xmlInputFileName, piroParams.get() );

      // =========================================================================
      // extract data of the parameter list
      Teuchos::ParameterList outputList = piroParams->sublist ( "Output", true );

      // set default directory to be the directory of the XML file itself
      std::string xmlPath = boost::filesystem::path ( xmlInputFileName ).branch_path().string();
      boost::filesystem::path outputDirectory = outputList.get<string> ( "Output directory" );
      if ( outputDirectory.root_directory().empty() ) // outputDirectory is empty or is a relative directory.
          outputDirectory = xmlPath / outputDirectory;
      std::string contFileBaseName =
          outputList.get<std::string> ( "Continuation file base name" );
      std::string outputFormat =
          outputList.get<std::string> ( "Output format" );
      boost::filesystem::path contDataFile =
          outputList.get<std::string> ( "Continuation data file name" );

      Teuchos::ParameterList initialGuessList;
      initialGuessList = piroParams->sublist ( "Initial guess", true );
      // =========================================================================

      Teuchos::ParameterList              problemParameters;
      Teuchos::RCP<Epetra_Vector>         z = Teuchos::null;
      Teuchos::RCP<VIO::EpetraMesh::Mesh> mesh = Teuchos::null;

      VIO::EpetraMesh::read( eComm,
                             getAbsolutePath( initialGuessList.get<std::string> ( "State" ), xmlPath ),
                             z,
                             mesh,
                             problemParameters
                           );

      // create the state
      TEUCHOS_ASSERT( !z.is_null() );
      Teuchos::RCP<Ginla::EpetraFVM::State> state =
              Teuchos::rcp( new Ginla::EpetraFVM::State( *z, mesh ) );

      // possibly overwrite the parameters
      Teuchos::ParameterList & overwriteParamsList = piroParams->sublist ( "Overwrite parameter list", true );
      bool overwriteParameters = overwriteParamsList.get<bool> ( "Overwrite parameters" );
      if ( overwriteParameters )
      {
          Teuchos::ParameterList & overwritePList = overwriteParamsList.sublist( "Parameters", true );
          problemParameters.setParameters( overwritePList );
      }

      double mu = problemParameters.get<double> ( "mu" );

      Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp =
              Teuchos::rcp ( new Ginla::MagneticVectorPotential::Z ( mu ) );

      // create the mode evaluator
      Teuchos::RCP<Ginla::EpetraFVM::ModelEvaluator> glModel =
                Teuchos::rcp( new Ginla::EpetraFVM::ModelEvaluator( mesh,
                                                                    problemParameters,
                                                                    mvp,
                                                                    state
                                                                  )
                            );

//      // check out args
//      std::cout << "FFF " << glModel->createOutArgs().supports( EpetraExt::ModelEvaluator::OUT_ARG_DfDp, 0 ).none() << std::endl;

      Teuchos::RCP<Ginla::IO::NoxObserver> observer;

      std::string contFilePath = getAbsolutePath( contDataFile, xmlPath );
      Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter =
          Teuchos::rcp( new Ginla::IO::StatsWriter( contFilePath ) );

      // warn if initial value was given twice
      std::string contParam = piroParams->sublist ( "LOCA" )
                                         .sublist ( "Stepper" )
                                         .get<std::string> ( "Continuation Parameter" );
      if ( problemParameters.isParameter ( contParam ) )
          std::cerr << "Warning: Continuation parameter \""
                    << contParam
                    << "\" explicitly given. Initial value will be overwritten by 'LOCA->Stepper->Initial Value', though."
                    << std::endl;

      // Use these two objects to construct a Piro solved application
      //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
      Teuchos::RCP<EpetraExt::ModelEvaluator> piro;

      // Declare the eigensaver; it will be used only for LOCA solvers, though.
      Teuchos::RCP<Ginla::IO::SaveEigenData> glEigenSaver;

      // Setup the data output.
      int maxLocaSteps = piroParams->sublist ( "LOCA" )
                                    .sublist ( "Stepper" )
                                    .get<int> ( "Max Steps" );
      Teuchos::RCP<Ginla::IO::StateWriter> stateWriter =
          Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory.string(),
                                                    "solution",
                                                    "vtu",
                                                    maxLocaSteps
                                                  )
                      );

      // switch by solver type
      std::string & solver = piroParams->get( "Piro Solver", "" );
      // ----------------------------------------------------------------------
      if ( solver == "NOX" )
      {
//          observer = Teuchos::rcp( new Ginla::IO::NoxObserver( stateWriter,
//                                                               glModel,
//                                                               Ginla::IO::NoxObserver::NONLINEAR,
//                                                               glModel
//                                                             )
//                                 );
//          observer->setStatisticsWriter( statsWriter );
//
//          piro = Teuchos::rcp(new Piro::Epetra::NOXSolver( piroParams,
//                                                           glModel,
//                                                           observer ));
      }
      // ----------------------------------------------------------------------
      else if ( solver == "LOCA" )
      {
          observer = Teuchos::null;
          observer = Teuchos::rcp( new Ginla::IO::NoxObserver( stateWriter,
                                                               glModel,
                                                               Ginla::IO::NoxObserver::CONTINUATION,
                                                               glModel
                                                             )
                                 );
          observer->setStatisticsWriter( statsWriter );

          Teuchos::RCP<LOCA::StatusTest::Combo> locaTest =
              Teuchos::rcp( new LOCA::StatusTest::Combo( LOCA::StatusTest::Combo::OR ) );

          // setup eingen saver
// #ifdef HAVE_LOCA_ANASAZI
//           Teuchos::ParameterList & eigenList = piroParams->sublist ( "LOCA" ).sublist ( "Stepper" ) .sublist ( "Eigensolver" );
//           std::string eigenvaluesFileName =
//               getAbsolutePath( outputList.get<std::string> ( "Eigenvalues file name" ), xmlPath );
//           std::string eigenstateFileNameAppendix =
//               outputList.get<std::string> ( "Eigenstate file name appendix" );
//
//           Teuchos::RCP<Ginla::IO::StatsWriter> eigenStatsWriter =
//               Teuchos::rcp( new Ginla::IO::StatsWriter( eigenvaluesFileName ) );
//
//           // initialize the stability change test with a pointer to the eigenvalue information
//           int stabilityChangeTreshold = 1; // stop when the stability changes by multiplicity 1
//           Teuchos::RCP<const Teuchos::ParameterList> eigendataList = eigenStatsWriter->getList();
//           Teuchos::RCP<LOCA::StatusTest::Abstract> stabilityChangeTest =
//               Teuchos::rcp( new Ginla::StatusTest::StabilityChange( eigendataList,
//                                                                     stabilityChangeTreshold ) );
//
//           locaTest->addStatusTest( stabilityChangeTest );
//
//           glEigenSaver = Teuchos::RCP<Ginla::IO::SaveEigenData> ( new Ginla::IO::SaveEigenData ( eigenList,
//                                                                                                  glModel,
//                                                                                                  stateWriter,
//                                                                                                  eigenStatsWriter ) );
//
//           Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy = glEigenSaver;
//           eigenList.set ( "Save Eigen Data Method", "User-Defined" );
//           eigenList.set ( "User-Defined Save Eigen Data Name", "glSaveEigenDataStrategy" );
//           eigenList.set ( "glSaveEigenDataStrategy", glSaveEigenDataStrategy );
// #endif

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
//                   Teuchos::RCP<VIO::Mesh::Mesh>    mesh = Teuchos::null;
//                   VIO::Mesh::read( Comm,
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
                                                          observer
//                                                           ,
//                                                           Teuchos::null,
//                                                           locaTest
                                                        ) );

//           // get stepper and inject it into the eigensaver
//           Teuchos::RCP<LOCA::Stepper> stepper = piroLOCASolver->getLOCAStepperNonConst();
// #ifdef HAVE_LOCA_ANASAZI
//           glEigenSaver->setLocaStepper ( stepper );
// #endif
          piro = piroLOCASolver;
      }
      // ----------------------------------------------------------------------
      else if ( solver == "Turning Point" )
      {
          // TODO make sure the turning point continuation doesn't technically fail by default

          observer = Teuchos::null;
//          observer = Teuchos::rcp( new Ginla::IO::NoxObserver( stateWriter,
//                                                               glModel,
//                                                               Ginla::IO::NoxObserver::TURNING_POINT ) );
//           observer->setStatisticsWriter( statsWriter, glOperator );

          // get the initial null state file
          initialGuessList = piroParams->sublist ( "Initial guess", true );

          // read the initial null state
          Teuchos::ParameterList              voidParameters;
          Teuchos::RCP<Epetra_Vector>         z = Teuchos::null;
          Teuchos::RCP<VIO::EpetraMesh::Mesh> mesh = Teuchos::null;
          VIO::EpetraMesh::read( eComm,
                                 getAbsolutePath( initialGuessList.get<string> ( "Null state" ), xmlPath),
                                 z,
                                 mesh,
                                 voidParameters
                               );

          TEUCHOS_ASSERT( !z.is_null() );
          Teuchos::RCP<Ginla::EpetraFVM::State> nullstate =
              Teuchos::rcp( new Ginla::EpetraFVM::State( *z, mesh ) );

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
          TEST_FOR_EXCEPTION( true,
                              std::logic_error,
                              "Unknown solver type \"" << solver << "\"." ) ;
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
      piro->evalModel(inArgs, outArgs);

      if ( outArgs.isFailed() )
          status += 10;

      // manually release LOCA stepper
#ifdef HAVE_LOCA_ANASAZI
      if ( !glEigenSaver.is_null() )
          glEigenSaver->releaseLocaStepper ();
#endif
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

    return status==0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
// =========================================================================
