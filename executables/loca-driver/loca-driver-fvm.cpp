#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

//#include <Piro_Epetra_NOXSolver.hpp>
//#include <Piro_Epetra_LOCASolver.hpp>
#include <LOCA_Epetra_Group.H>
#include <LOCA_Thyra_Group.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <LOCA_Epetra_ModelEvaluatorInterface.H>
#include <LOCA_Epetra_Factory.H>
#include <NOX_StatusTest_Generic.H>
#include <NOX_StatusTest_Factory.H>
#include <LOCA_StatusTest_Combo.H>
#include <LOCA_StatusTest_MaxIters.H>

#include <boost/filesystem.hpp>

#include "Ginla_EpetraFVM_ModelEvaluator.h"
#include "Ginla_EpetraFVM_State.h"
#include "VIO_EpetraMesh_Reader.h"

// #include "Ginla_IO_SaveNewtonData.h"
#include "Ginla_IO_SaveEigenData.h"
#include "Ginla_IO_NoxObserver.h"

#include "Ginla_Helpers.h"

#include "Ginla_MagneticVectorPotential_X.h"
#include "Ginla_MagneticVectorPotential_Y.h"
#include "Ginla_MagneticVectorPotential_Z.h"

#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_StatsWriter.h"

#include "Ginla_StatusTest_MaxAcceptedSteps.h"
#include "Ginla_StatusTest_Energy.h"
#include "Ginla_StatusTest_Loop.h"
#include "Ginla_StatusTest_ParameterLimits.h"
#include "Ginla_StatusTest_StabilityChange.h"



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
int
main ( int argc, char *argv[] )
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
        // =====================================================================
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

        // =====================================================================
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
        // =====================================================================

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


        // set I/O routines
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
        std::string contFilePath = getAbsolutePath( contDataFile, xmlPath );
        Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter =
                Teuchos::rcp( new Ginla::IO::StatsWriter( contFilePath ) );

        Teuchos::RCP<Ginla::IO::NoxObserver> observer =
                Teuchos::rcp( new Ginla::IO::NoxObserver( stateWriter,
                                                          glModel,
                                                          Ginla::IO::NoxObserver::CONTINUATION,
                                                          glModel
                                                        )
                            );
        observer->setStatisticsWriter( statsWriter );


        // explicitly copy the parameter list as LOCA::createGlobalData needs a non-const
        Teuchos::RCP<Teuchos::ParameterList> locaParams =
                Teuchos::rcp( new Teuchos::ParameterList( piroParams ->sublist ( "LOCA" ) ) );

        std::string contParam = piroParams->sublist ( "LOCA" )
                                .sublist ( "Stepper" )
                                .get<std::string> ( "Continuation Parameter" );
        if ( problemParameters.isParameter ( contParam ) )
            std::cerr << "Warning: Continuation parameter \""
                    << contParam
                    << "\" explicitly given. Initial value will be overwritten by 'LOCA->Stepper->Initial Value', though."
                    << std::endl;

        // ---------------------------------------------------------------------
        // Create Epetra factory
        Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
                Teuchos::rcp ( new LOCA::Epetra::Factory() );

        // Create global data object
        Teuchos::RCP<LOCA::GlobalData> globalData =
                LOCA::createGlobalData ( locaParams, epetraFactory );
        // ---------------------------------------------------------------------
        Teuchos::RCP<LOCA::Epetra::ModelEvaluatorInterface> locaModelEvaluatorInterface
                = Teuchos::rcp( new LOCA::Epetra::ModelEvaluatorInterface( globalData, glModel ) );
        // ---------------------------------------------------------------------
        locaModelEvaluatorInterface->setObserver( observer );
        // ---------------------------------------------------------------------
        // Create all possible Epetra_Operators.
        Teuchos::RCP<Epetra_Operator> J = glModel->create_W();
        Teuchos::RCP<Epetra_Operator> M = glModel->create_WPrec()->PrecOp;


        Epetra_Vector X( z->Map() );
        X.Random();
        Epetra_Vector Y( X );
        M->Apply( X, Y );
        std::cout << "M successfully applied" << std::endl;
        M->ApplyInverse( X, Y );
        std::cout << "M successfully ApplyInversed" << std::endl;

        // Create the linear system.
        // Use the TimeDependent interface for computation of shifted matrices.
        //Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = locaModelEvaluatorInterface;
        Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = locaModelEvaluatorInterface;
        Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = locaModelEvaluatorInterface;

        Teuchos::ParameterList& nlPrintParams =
                piroParams->sublist ( "NOX" ) .sublist ( "Printing", true );

        Teuchos::ParameterList& lsParams =
                piroParams->sublist ( "NOX" ) .sublist ( "Direction" ) .sublist ( "Newton" ) .sublist ( "Linear Solver", true );

        //  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(
        //      new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams, iJac, J, iPrec, M, *soln));

        NOX::Epetra::Vector cloneVector( Epetra_Vector( *glModel->get_x_map() ) );
        Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
                Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( nlPrintParams,
                                                                      lsParams,
                                                                      iJac,
                                                                      J,
                                                                      iPrec,
                                                                      M,
                                                                      cloneVector ) );

        Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = locaModelEvaluatorInterface;
        // ---------------------------------------------------------------------
        // Create a group which uses that problem interface. The group will
        // be initialized to contain the default initial guess for the
        // specified problem.
        LOCA::ParameterVector myLocaParams = *(glModel->getParameters());

        // get initial guess
        NOX::Epetra::Vector initialGuess ( state->getPsiNonConst(),
                                           NOX::Epetra::Vector::CreateView
                                           );

        Teuchos::RCP<LOCA::Epetra::Group> grp =
                Teuchos::rcp ( new LOCA::Epetra::Group ( globalData,
                                                         nlPrintParams,
                                                         iTime,
                                                         initialGuess,
                                                         linSys,
                                                         linSys,
                                                         myLocaParams ) );

        std::cout << myLocaParams << std::endl;
        grp->setParams ( myLocaParams );
        // ---------------------------------------------------------------------
        // Set up the NOX status tests
        Teuchos::ParameterList& noxList = piroParams->sublist ( "NOX", true );
        Teuchos::RCP<NOX::StatusTest::Generic> noxStatusTest = NOX::StatusTest::buildStatusTests( noxList.sublist( "Status Tests" ),
                                                                                                  *globalData->locaUtils()
                                                                                                  );
        // ---------------------------------------------------------------------
        // Set up the LOCA status tests
        Teuchos::RCP<LOCA::StatusTest::Combo> locaCombo =
                Teuchos::rcp ( new LOCA::StatusTest::Combo ( LOCA::StatusTest::Combo::OR ) );

        Teuchos::RCP<LOCA::StatusTest::Abstract> maxLocaStepsTest =
                Teuchos::rcp ( new LOCA::StatusTest::MaxIters ( maxLocaSteps ) );
        locaCombo->addStatusTest( maxLocaStepsTest );

        //      Teuchos::RCP<LOCA::StatusTest::Abstract> zeroEnergyTest =
        //                       Teuchos::rcp ( new Ginla::StatusTest::Energy ( glModel,
        //                                                                      0.0 ) );

        //      Teuchos::RCP<LOCA::StatusTest::Abstract> loopTest =
        //                       Teuchos::rcp ( new Ginla::StatusTest::Loop ( glsystem ) );

        //      Teuchos::RCP<LOCA::StatusTest::Abstract> turnaroundTest =
        //                       Teuchos::rcp ( new Ginla::StatusTest::Turnaround () );



        //     locaCombo->addStatusTest( zeroEnergyTest );
        //      locaCombo->addStatusTest( loopTest );
        //      locaCombo->addStatusTest( turnaroundTest );
        // ---------------------------------------------------------------------
        // Create the stepper
        Teuchos::RCP<LOCA::Stepper> stepper =
                Teuchos::rcp ( new LOCA::Stepper ( globalData,
                                                   grp,
                                                   locaCombo,
                                                   noxStatusTest,
                                                   piroParams ) );
        // ---------------------------------------------------------------------
        // pass pointer to stepper to glsystem to be able to read stats from the stepper in there
//        locaModelEvaluatorInterface->setLocaStepper ( stepper );
        // ---------------------------------------------------------------------
        // Perform continuation run
        status = stepper->run();
        // ---------------------------------------------------------------------
        // clean up
        LOCA::destroyGlobalData ( globalData );
//        locaModelEvaluatorInterface->releaseLocaStepper();
        // ---------------------------------------------------------------------
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