// see <http://old.nabble.com/Undefined-reference-to-%27main%27-with-Boost-Test.-Why--td15986217.html>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GRNN zeroStepLocaTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Epetra_CrsMatrix.h>

#include <LOCA_Epetra_Factory.H>
#include <LOCA_Epetra_Group.H>
#include <LOCA_StatusTest_MaxIters.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_StatusTest_NormF.H>
#include <NOX_StatusTest_MaxIters.H>
#include <NOX_StatusTest_Combo.H>
// #include <NOX_Epetra_Group.H>

#include "Ginla_MagneticVectorPotential_ZSquareSymmetric.h"
#include "Ginla_FDM_Operator_BCCentral.h"
#include "Ginla_IO_StatsWriter.h"
#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_SaveEigenData.h"
#include "Ginla_Helpers.h"
#include "Ginla_FDM_LocaSystem_Bordered.h"

#include "Recti_Grid_Uniform.h"
#include "Recti_Grid_Reader.h"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <boost/filesystem.hpp>
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// ===========================================================================
// <http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/tutorials/hello-the-testing-world.htmlhh>
BOOST_AUTO_TEST_CASE( zero_step_loca_test )
{
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( &boost::unit_test::framework::master_test_suite().argc,
               &boost::unit_test::framework::master_test_suite().argv
             );
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

    // ------------------------------------------------------------------------
    // handle command line arguments
    Teuchos::CommandLineProcessor My_CLP;

    std::string xmlInputFileName = "";
    My_CLP.setOption ( "xml-input-file", &xmlInputFileName,
                       "XML file containing the parameter list", true );
                       
    std::string expSolFileName = "";
    My_CLP.setOption ( "expected-solution-file", &expSolFileName,
                       "VTK/VTI file containing the expected solution", true );

    // print warning for unrecognized arguments
    My_CLP.recogniseAllOptions ( true );

    // don't throw exceptions
    My_CLP.throwExceptions ( true );

    // finally, parse the stuff!
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
    
    parseReturn = My_CLP.parse ( boost::unit_test::framework::master_test_suite().argc,
                                 boost::unit_test::framework::master_test_suite().argv );
    // ------------------------------------------------------------------------
    // Read the XML file
    Teuchos::RCP<Teuchos::ParameterList> paramList =
            Teuchos::rcp ( new Teuchos::ParameterList );

    Teuchos::updateParametersFromXmlFile ( xmlInputFileName,
                                           paramList.get()
                                         );
    // ------------------------------------------------------------------------
    // extract data for the output the parameter list
    Teuchos::ParameterList outputList;
    outputList = paramList->sublist ( "Output", true );

    // set default directory to be the directory of the XML file itself
    std::string xmlPath = boost::filesystem::path ( xmlInputFileName ).branch_path().string();
    boost::filesystem::path outputDirectory = outputList.get<string> ( "Output directory" );
    if ( outputDirectory.root_directory().empty() ) // outputDirectory is empty or is a relative directory.
        outputDirectory = xmlPath / outputDirectory;
    std::string contFileBaseName =
        outputList.get<string> ( "Continuation file base name" );
    std::string outputFormat =
        outputList.get<string> ( "Output format" );
    std::string contDataFileName =
        outputList.get<string> ( "Continuation data file name" );
    // ------------------------------------------------------------------------
    Teuchos::ParameterList initialGuessList;
    initialGuessList = paramList->sublist ( "Initial guess", true );

    boost::filesystem::path inputGuessFile = initialGuessList.get<string> ( "File name" );
    if ( !inputGuessFile.empty() && inputGuessFile.root_directory().empty() ) // if inputGuessFile is a relative path
        inputGuessFile = xmlPath / inputGuessFile;

    BOOST_REQUIRE( !inputGuessFile.empty() );
    
    // For technical reasons, the reader can only accept ComplexMultiVectors.
    Teuchos::RCP<Ginla::FDM::State> state;
    Teuchos::RCP<Recti::Grid::Uniform> grid;
    Teuchos::ParameterList glParameters;
    
    Recti::Grid::Reader::read ( Comm, inputGuessFile.string(), state, grid, glParameters );

    // possibly overwrite the parameters
    Teuchos::ParameterList & overwriteParamsList = paramList->sublist ( "Overwrite parameter list", true ); 
    bool overwriteParameters = overwriteParamsList.get<bool> ( "Overwrite parameters" );
    if ( overwriteParameters )
    {
        Teuchos::ParameterList & overwritePList = overwriteParamsList.sublist( "Parameters", true );
        glParameters.setParameters( overwritePList );
        
        // possibly update the scaling of the grid
        grid->updateScaling( glParameters.get<double>("scaling") );
    }
    // ------------------------------------------------------------------------

    Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::ZSquareSymmetric
                                  ( glParameters.get<double> ( "H0" ),
                                    glParameters.get<double> ( "scaling" )
                                  )
                     );

    // create the operator
    Teuchos::RCP<const ComplexMap> map = state->getPsi()->getMap();
    Teuchos::RCP<Ginla::FDM::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::FDM::Operator::BCCentral ( grid, A, map, map ) );

    std::string fn = outputDirectory.string() + "/" + contDataFileName;
    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter = 
            Teuchos::rcp( new Ginla::IO::StatsWriter( fn ) );

    Teuchos::ParameterList & stepperList = paramList->sublist ( "LOCA" ).sublist ( "Stepper" );
    int maxLocaSteps = stepperList.get<int> ( "Max Steps" );
        
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter =
        Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory.string(),
                                                  contFileBaseName,
                                                  outputFormat,
                                                  maxLocaSteps ) );
    
    Teuchos::RCP<Ginla::FDM::LocaSystem::Bordered> glsystem =
               Teuchos::rcp ( new Ginla::FDM::LocaSystem::Bordered ( glOperator,
                                                                eComm,
                                                                map,
                                                                statsWriter,
                                                                stateWriter ) );
    
    // set the initial value from glParameters
    std::string contParam = stepperList.get<string> ( "Continuation Parameter" );

    BOOST_REQUIRE( glParameters.isParameter ( contParam ) );

    // check if the initial value was given (will be unused anyway)
    if ( stepperList.isParameter ( "Initial Value" ) )
    {
        std::cerr << "Warning: Parameter 'LOCA->Stepper->Initial Value' given, but will not be used."
                  << std::endl;
    }

    stepperList.set ( "Initial Value", glParameters.get<double> ( contParam ) );

    // ------------------------------------------------------------------------
    // Create the necessary objects
    // ------------------------------------------------------------------------
    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
        Teuchos::rcp ( new LOCA::Epetra::Factory );

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
        LOCA::createGlobalData ( paramList, epetraFactory );
    // ------------------------------------------------------------------------
#ifdef HAVE_LOCA_ANASAZI
    Teuchos::ParameterList& eigenList =
        paramList->sublist ( "LOCA" ).sublist ( "Stepper" ) .sublist ( "Eigensolver" );
    std::string eigenvaluesFileName =
        outputList.get<string> ( "Eigenvalues file name" );
    std::string eigenstateFileNameAppendix =
        outputList.get<string> ( "Eigenstate file name appendix" );

    Teuchos::RCP<Ginla::IO::StatsWriter> eigenStatsWriter =
        Teuchos::rcp( new Ginla::IO::StatsWriter( eigenvaluesFileName ) );

    Teuchos::RCP<Ginla::IO::SaveEigenData> glEigenSaver =    
        Teuchos::RCP<Ginla::IO::SaveEigenData> ( new Ginla::IO::SaveEigenData ( eigenList,
                                                                                glsystem,
                                                                                stateWriter,
                                                                                eigenStatsWriter ) );

    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy = glEigenSaver;
    eigenList.set ( "glSaveEigenDataStrategy", glSaveEigenDataStrategy );
#endif

    // ------------------------------------------------------------------------
    // Create all possible Epetra_Operators.
    Teuchos::RCP<Epetra_RowMatrix> J = glsystem->getJacobian();
//  Teuchos::RCP<Epetra_RowMatrix> M = glsystem->getPreconditioner();

    // Create the linear system.
    // Use the TimeDependent interface for computation of shifted matrices.
    Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = glsystem;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glsystem;
//  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = glsystem;

    Teuchos::ParameterList& nlPrintParams =
        paramList->sublist ( "NOX" )
                  .sublist ( "Printing" );

    Teuchos::ParameterList& lsParams =
        paramList->sublist ( "NOX" )
                  .sublist ( "Direction" )
                  .sublist ( "Newton" )
                  .sublist ( "Linear Solver" );

    NOX::Epetra::Vector cloneVector( Epetra_Vector( *glsystem->getMap() ) );
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
        Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( nlPrintParams,
                                                              lsParams,
                                                              iReq,
                                                              iJac,
                                                              J,
                                                              cloneVector ) );

    Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = glsystem;
    // ------------------------------------------------------------------------
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.

    // translate parameters into a LOCA list
    LOCA::ParameterVector locaParams =
          *(Ginla::Helpers::teuchosParameterList2locaParameterVector( glParameters ));

    // get initial guess
    NOX::Epetra::Vector initialGuess ( glsystem->createSystemVector( *state ),
                                       NOX::Epetra::Vector::CreateView
                                     );

    Teuchos::RCP<LOCA::Epetra::Group> grp =
        Teuchos::rcp ( new LOCA::Epetra::Group ( globalData,
                                                 nlPrintParams,
                                                 iTime,
                                                 initialGuess,
                                                 linSys,
                                                 linSys,
                                                 locaParams ) );

    grp->setParams ( locaParams );
    // ------------------------------------------------------------------------
    // Set up the NOX status tests
    Teuchos::ParameterList& noxList = paramList->sublist ( "NOX", true );
    double tol = noxList.get<double> ( "Tolerance" );
    int maxNonlinarSteps = noxList.get<int> ( "Max steps" );
    Teuchos::RCP<NOX::StatusTest::NormF> normF =
        Teuchos::rcp ( new NOX::StatusTest::NormF ( tol ) );
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters =
        Teuchos::rcp ( new NOX::StatusTest::MaxIters ( maxNonlinarSteps ) );
    Teuchos::RCP<NOX::StatusTest::Generic> comboOR =
        Teuchos::rcp ( new NOX::StatusTest::Combo ( NOX::StatusTest::Combo::OR,
                                                    normF,
                                                    maxIters ) );
    // ------------------------------------------------------------------------
    // Set up the LOCA status tests
    Teuchos::RCP<LOCA::StatusTest::MaxIters> maxLocaStepsTest;
    
    maxLocaStepsTest = Teuchos::rcp ( new LOCA::StatusTest::MaxIters ( maxLocaSteps ) );
    // ------------------------------------------------------------------------
    // Create the stepper
    Teuchos::RCP<LOCA::Stepper> stepper =
        Teuchos::rcp ( new LOCA::Stepper ( globalData, grp, maxLocaStepsTest, comboOR, paramList ) );
    // ------------------------------------------------------------------------
    // make sure that the stepper starts off with the correct starting value
    // pass pointer to stepper to glsystem to be able to read stats from the stepper in there
    glsystem->setLocaStepper ( stepper );
#ifdef HAVE_LOCA_ANASAZI
    glEigenSaver->setLocaStepper ( stepper );
#endif
    // ------------------------------------------------------------------------
    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status;
    status = stepper->run();
    // ------------------------------------------------------------------------
    // retrieve solution
    const NOX::Epetra::Group & finalGroup =
        Teuchos::dyn_cast<const NOX::Epetra::Group> ( *stepper->getSolutionGroup() );

    const Epetra_Vector & finalSolution =
        ( Teuchos::dyn_cast<const NOX::Epetra::Vector> ( finalGroup.getX() ) ).getEpetraVector();
    Teuchos::RCP<Ginla::State::Virtual> solutionState = glsystem->createState( finalSolution );
    // ------------------------------------------------------------------------
    // read expected solution from file
    // For technical reasons, the reader can only accept ComplexMultiVectors.
    Teuchos::RCP<Ginla::FDM::State> referenceState;
    Recti::Grid::Reader::read ( Comm, expSolFileName, referenceState, grid, glParameters );
    // ------------------------------------------------------------------------
    // compare the results:
    // get final solution
    solutionState->update( -1.0, *referenceState, 1.0 );
    
    double nrm = solutionState->normalizedScaledL2Norm();

    BOOST_CHECK_SMALL( nrm, 1.0e-10 );
    // ------------------------------------------------------------------------
    // clean up
    LOCA::destroyGlobalData ( globalData );
    glsystem->releaseLocaStepper();
#ifdef HAVE_LOCA_ANASAZI
    glEigenSaver->releaseLocaStepper();
#endif
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    // ------------------------------------------------------------------------

    return;
}
// ============================================================================
