// see <http://old.nabble.com/Undefined-reference-to-%27main%27-with-Boost-Test.-Why--td15986217.html>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GRNN zeroStepPiroTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

#include "Recti_Grid_Uniform.h"
#include "Recti_Grid_Reader.h"

// #include "Ginla_IO_SaveNewtonData.h"
#include "Ginla_IO_SaveEigenData.h"
#include "Ginla_IO_NoxObserver.h"
#include "Ginla_Komplex_LinearProblem.h"

#include "Ginla_FDM_ModelEvaluator_Default.h"
#include "Ginla_FDM_ModelEvaluator_Bordered.h"
#include "Ginla_FDM_Operator_BCCentral.h"
// #include "Ginla_Operator_BCInner.h"
// #include "Ginla_Operator_BCOuter.h"

#include "Ginla_MagneticVectorPotential_Centered.h"
#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_StatsWriter.h"
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// ===========================================================================
// <http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/tutorials/hello-the-testing-world.htmlhh>
BOOST_AUTO_TEST_CASE( branch_switch_test )
{
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init( boost::unit_test::framework::master_test_suite().argc,
              boost::unit_test::framework::master_test_suite().argv );
#endif

    // Create a communicator for Tpetra objects
    const Teuchos::RCP<const Teuchos::Comm<int> > Comm =
          Teuchos::DefaultComm<int>::getComm();

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm =
      Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm>  eComm =
           Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

    // ===========================================================================
    // handle command line arguments
    Teuchos::CommandLineProcessor My_CLP;

    My_CLP.setDocString (
        "This program solves the Ginzburg--Landau problem with a NOX interface.\n"
    );

    std::string xmlInputFileName = "";
    My_CLP.setOption ( "xml-input-file", &xmlInputFileName,
                      "XML file containing the parameter list", true );

    std::string expSolFileName = "";
    My_CLP.setOption ( "expected-solution-file", &expSolFileName,
                      "VTK/VTI file containing the expected solution", true );
                      
    // print warning for unrecognized arguments
    My_CLP.recogniseAllOptions ( true );

    // do throw exceptions
    My_CLP.throwExceptions ( true );

    int status = 0;
    
    // finally, parse the command line
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
    parseReturn = My_CLP.parse ( boost::unit_test::framework::master_test_suite().argc,
                                 boost::unit_test::framework::master_test_suite().argv );

    Teuchos::RCP<Teuchos::ParameterList> piroParams =
        Teuchos::rcp ( new Teuchos::ParameterList );
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

    Teuchos::ParameterList             problemParameters;
    Teuchos::RCP<Ginla::FDM::State>         state;
    Teuchos::RCP<Recti::Grid::Uniform> grid = Teuchos::null;

    boost::filesystem::path statefile = initialGuessList.get<std::string> ( "State" );
    if ( !statefile.empty() && statefile.root_directory().empty() ) // filepath is a relative path
        statefile = xmlPath / statefile;
    TEUCHOS_ASSERT( !statefile.empty() );

    Recti::Grid::Reader::read ( Comm,
                                statefile.string(),
                                state,
                                grid,
                                problemParameters );

    // possibly overwrite the parameters
    Teuchos::ParameterList & overwriteParamsList = piroParams->sublist ( "Overwrite parameter list", true ); 
    bool overwriteParameters = overwriteParamsList.get<bool> ( "Overwrite parameters" );
    if ( overwriteParameters )
    {
        Teuchos::ParameterList & overwritePList = overwriteParamsList.sublist( "Parameters", true );
        problemParameters.setParameters( overwritePList );
        
        // possibly update the scaling of the grid
        grid->updateScaling( problemParameters.get<double>("scaling") );
    }
                                
    double h0      = problemParameters.get<double> ( "H0" );
    double scaling = problemParameters.get<double> ( "scaling" );
    
    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( h0 ) );

    Teuchos::ParameterList & stepperList = piroParams->sublist ( "LOCA" ).sublist ( "Stepper" );
    int maxLocaSteps = stepperList.get<int> ( "Max Steps" );
        
    // setup the data output
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter =
        Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory.string(),
                                                  "solution",
                                                  "VTI",
                                                  maxLocaSteps ) );

    // create the operator
    Teuchos::RCP<Ginla::FDM::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::FDM::Operator::BCCentral ( grid,
                                                        A,
                                                        state->getPsi()->getMap(),
                                                        state->getPsi()->getMap() ) );
                                                        
    Teuchos::RCP<Ginla::Komplex::LinearProblem> komplex =
        Teuchos::rcp( new Ginla::Komplex::LinearProblem( eComm, state->getPsi()->getMap() ) );

    // create the mode evaluator
    Teuchos::RCP<Ginla::FDM::ModelEvaluator::Default> glModel = 
              Teuchos::rcp(new Ginla::FDM::ModelEvaluator::Default( glOperator,
                                                               komplex,
                                                               *state,
                                                               problemParameters ) );
                                                              
    Teuchos::RCP<Ginla::IO::NoxObserver> observer;
                                                  
    std::string contFilePath = (outputDirectory / contDataFile).string();
    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter =
        Teuchos::rcp( new Ginla::IO::StatsWriter( contFilePath ) );
    
    // set the initial value from glParameters
    std::string contParam = stepperList.get<std::string> ( "Continuation Parameter" );
    TEST_FOR_EXCEPTION ( !problemParameters.isParameter ( contParam ),
                        std::logic_error,
                        "Parameter \"" << contParam << "\" given as continuation parameter, but doesn't exist"
                        << "in the glParameters list." );
                        
    // check if the initial value was given (won't be unused anyway)
    if ( stepperList.isParameter ( "Initial Value" ) )
        std::cerr << "Warning: Parameter 'LOCA->Stepper->Initial Value' given, but will not be used."
                  << std::endl;

    // Set initial value appropriately.
    // TODO Get rid of the explicit "double".
    stepperList.set ( "Initial Value", problemParameters.get<double> ( contParam ) );

    
    // Use these two objects to construct a Piro solved application 
    //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
    Teuchos::RCP<EpetraExt::ModelEvaluator> piro;
    
    // Declare the eigensaver; it will be used only for LOCA solvers, though.
    Teuchos::RCP<Ginla::IO::SaveEigenData> glEigenSaver;
    
    // switch by solver type
    std::string & solver = piroParams->get( "Piro Solver", "" );


    observer = Teuchos::rcp( new Ginla::IO::NoxObserver( stateWriter,
                                                         glModel,
                                                         Ginla::IO::NoxObserver::CONTINUATION ) );
    observer->setStatisticsWriter( statsWriter, glOperator );

    // feed in restart predictor vector
    Teuchos::ParameterList & predictorList = piroParams->sublist ( "LOCA" )
                                                        .sublist ( "Predictor" );
    if ( predictorList.get<std::string>("Method") == "Secant" )
    {
        Teuchos::ParameterList & fspList = predictorList.sublist ( "First Step Predictor" );
        if ( fspList.get<std::string>("Method") == "Restart" )
        { 
            boost::filesystem::path predictorfile = initialGuessList.get<std::string> ( "Predictor" );
            if ( !predictorfile.empty() && predictorfile.root_directory().empty() ) // filepath is a relative path
                predictorfile = xmlPath / predictorfile;
            TEUCHOS_ASSERT( !predictorfile.empty() );
          
            // Get restart vector.
            // read the predictor state
            Teuchos::ParameterList             voidParameters;
            Teuchos::RCP<Ginla::FDM::State>         predictorState;
            Teuchos::RCP<Recti::Grid::Uniform> grid = Teuchos::null;
            Recti::Grid::Reader::read ( Comm,
                                        predictorfile.string(),
                                        predictorState,
                                        grid,
                                        voidParameters );
                                        
            // transform to system vector
            Teuchos::RCP<Epetra_Vector> predictorV = glModel->createSystemVector( *predictorState );
                                                          
            // Create the LOCA::MultiContinuation::ExtendedVector;
            // that's the predictor vector with a prediction for the parameter appended
            // -- 0 in the case of a pitchfork.
            NOX::Epetra::Vector vx( predictorV );
            
            // Note that the globalData element is here the null pointer. It would be work-aroundish to
            // artificially create one here, and after all, all that's it's used for is error and warning
            // message printing.
            Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> restartVector =
                Teuchos::rcp( new LOCA::MultiContinuation::ExtendedVector( Teuchos::null,
                                                                            vx,
                                                                            1 )
                            );
            double vp = 0.0;
            restartVector->setScalar(0, vp);
            
            fspList.set( "Restart Vector", restartVector );
        }
    }

    // fetch the stepper
    piro =  Teuchos::rcp(new Piro::Epetra::LOCASolver( piroParams,
                                                       glModel,
                                                       observer ));


    // Now the (somewhat cumbersome) setting of inputs and outputs
    EpetraExt::ModelEvaluator::InArgs inArgs = piro->createInArgs();
    int num_p = inArgs.Np();     // Number of *vectors* of parameters
    Teuchos::RCP<Epetra_Vector> p1 =
        Teuchos::rcp(new Epetra_Vector(*(piro->get_p_init(0))));
    inArgs.set_p(0,p1);   

    // Set output arguments to evalModel call
    EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();

    // Now, solve the problem and return the responses
    piro->evalModel(inArgs, outArgs);
    
    // ------------------------------------------------------------------------
    // fetch the solution
    // outArgs.get_g(0) must be gx
    BOOST_ASSERT( !outArgs.get_g(0).is_null() ); 
    Teuchos::RCP<Ginla::State::Virtual> solutionState = glModel->createState( *(outArgs.get_g(0)) );
    // ------------------------------------------------------------------------
    // read reference solution from file
    Teuchos::RCP<Ginla::FDM::State> refState;
    Recti::Grid::Reader::read ( Comm, expSolFileName, refState, grid, problemParameters );
    // ------------------------------------------------------------------------
    // compare the results:
    // get final solution
    
//     solutionState->save( "test.vti" );
    
    Teuchos::RCP<Ginla::State::Virtual> diff = solutionState;
    
    diff->update( -1.0, *refState, 1.0 );    
    
    // don't be as strict here: 10^{-5} is enough
    BOOST_CHECK_SMALL( diff->normalizedScaledL2Norm(), 1.0e-5 );
    // ------------------------------------------------------------------------

#ifdef HAVE_MPI
      MPI_Finalize();
#endif

    return;
}
// ============================================================================
