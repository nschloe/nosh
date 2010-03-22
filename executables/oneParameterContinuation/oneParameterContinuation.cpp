#include <boost/filesystem.hpp>

#include <LOCA.H>
#include <LOCA_Epetra_Factory.H>
#include <LOCA_Epetra_Group.H>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Vector.h>
#include <NOX_Epetra_LinearSystem_AztecOO.H>

#include <LOCA_StatusTest_MaxIters.H>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "Ginla_Operator_BCCentral.h"
#include "Ginla_Operator_BCInner.h"

#include "Recti_Domain_Square.h"

#include "Ginla_LocaSystem_Bordered.h"
#include "Ginla_IO_SaveEigenData.h"
#include "Ginla_Helpers.h"
#include "Ginla_IO_StatsWriter.h"
#include "Ginla_IO_StateWriter.h"

#include "Recti_Grid_Reader.h"

#include "Ginla_Perturbation_Quadrants.h"

#include <LOCA_Thyra_Group.H>

// =============================================================================
int
main ( int argc, char *argv[] )
{
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( &argc,&argv );
#endif

    // create Epetra communicator
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm =
        Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm> eComm =
        Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

    // Create a communicator for Tpetra objects
    const Teuchos::RCP<const Teuchos::Comm<int> > Comm =
        Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // handle command line arguments
    Teuchos::CommandLineProcessor My_CLP;

    My_CLP.setDocString (
        "This program does continuation for the Ginzburg--Landau problem with a LOCA interface.\n"
        "It is possible to give an initial guess in VTK format on the command line.\n" );

    std::string xmlInputFileName = "";
    My_CLP.setOption ( "xml-input-file", &xmlInputFileName,
                       "XML file containing the parameter list", true );

    // print warning for unrecognized arguments
    My_CLP.recogniseAllOptions ( true );

    // don't throw exceptions
    My_CLP.throwExceptions ( false );

    // finally, parse the stuff!
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
    try
    {
        parseReturn = My_CLP.parse ( argc, argv );
    }
    catch ( std::exception &e )
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    Teuchos::RCP<Teuchos::ParameterList> paramList =
        Teuchos::rcp ( new Teuchos::ParameterList );
    std::cout << "Reading parameter list from \"" << xmlInputFileName << "\"."
              << std::endl;

    try
    {
        Teuchos::updateParametersFromXmlFile ( xmlInputFileName, paramList.get() );
    }
    catch ( std::exception &e )
    {
        std::cerr << e.what() << "\n"
                  << "Check your XML file for syntax errors." << std::endl;
        return 1;
    }
    // =========================================================================
    // extract data for the output the parameter list
    Teuchos::ParameterList outputList;
    try
    {
        outputList = paramList->sublist ( "Output", true );
    }
    catch ( std::exception &e )
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
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
    // =========================================================================


    // ---------------------------------------------------------------------------
    Teuchos::ParameterList glParameters;
    Teuchos::RCP<ComplexVector> psi;
    Teuchos::RCP<Recti::Grid::Uniform> grid;

    Teuchos::ParameterList initialGuessList;
    try
    {
        initialGuessList = paramList->sublist ( "Initial guess", true );
    }
    catch ( std::exception &e )
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    boost::filesystem::path inputGuessFile = initialGuessList.get<string> ( "File name" );
    if ( !inputGuessFile.empty() && inputGuessFile.root_directory().empty() ) // if inputGuessFile is a relative path
        inputGuessFile = xmlPath / inputGuessFile;

    TEUCHOS_ASSERT( !inputGuessFile.empty() );
    
    try
    {
        // For technical reasons, the reader can only accept ComplexMultiVectors.
        Teuchos::RCP<ComplexMultiVector> psiM;
        Recti::Grid::Reader::read ( Comm, inputGuessFile.string(), psiM, grid, glParameters );
        TEUCHOS_ASSERT_EQUALITY ( psiM->getNumVectors(), 1 );
        psi = psiM->getVectorNonConst ( 0 );
    }
    catch ( std::exception &e )
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }

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
    // ---------------------------------------------------------------------------

    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( glParameters.get<double> ( "H0" ),
                                                                      glParameters.get<double> ( "scaling" ) ) );

    // create the operator
    Teuchos::RCP<Ginla::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::Operator::BCCentral ( grid, A ) );

//     // create a perturbation
//     Teuchos::RCP<GL::Perturbation::Virtual> quadrantsPerturbation =
//       Teuchos::rcp ( new GL::Perturbation::Quadrants ( grid ) );

    // TODO: why not make glOperator depend upon perturbation instead?
//     GinzburgLandau glProblem = GinzburgLandau ( glOperator,
//                                                 quadrantsPerturbation );


    std::string fn = outputDirectory.string() + "/" + contDataFileName;
    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter = 
        Teuchos::rcp( new Ginla::IO::StatsWriter( fn ) );

    Teuchos::RCP<Ginla::LocaSystem::Bordered> glsystem;

    Teuchos::ParameterList & stepperList = paramList->sublist ( "LOCA" ).sublist ( "Stepper" );
    int maxLocaSteps = stepperList.get<int> ( "Max Steps" );
    
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter = 
        Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory.string(),
                                                  contFileBaseName,
                                                  outputFormat,
                                                  maxLocaSteps ) );
    try
    {
        glsystem = Teuchos::rcp ( new Ginla::LocaSystem::Bordered ( glOperator,
                                                                    eComm,
                                                                    psi->getMap(),
                                                                    statsWriter,
                                                                    stateWriter ) );
    }
    catch ( std::exception & e )
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    
    // set the initial value from glParameters
    std::string contParam = stepperList.get<string> ( "Continuation Parameter" );
    TEST_FOR_EXCEPTION ( !glParameters.isParameter ( contParam ),
                         std::logic_error,
                         "Parameter \"" << contParam << "\" given as continuation parameter, but doesn't exist"
                         << "in the glParameters list." );

    // check if the initial value was given (will be unused anyway)
    if ( stepperList.isParameter ( "Initial Value" ) )
    {
        std::cerr << "Warning: Parameter 'LOCA->Stepper->Initial Value' given, but will not be used."
                  << std::endl;
    }

    // TODO Get rid of the explicit "double".
    stepperList.set ( "Initial Value", glParameters.get<double> ( contParam ) );

    // ---------------------------------------------------------------------------
    // Create the necessary objects
    // ---------------------------------------------------------------------------
    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
        Teuchos::rcp ( new LOCA::Epetra::Factory );

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
        LOCA::createGlobalData ( paramList, epetraFactory );
    // ---------------------------------------------------------------------------

#ifdef HAVE_LOCA_ANASAZI
    Teuchos::ParameterList& eigenList =
        paramList->sublist ( "LOCA" ).sublist ( "Stepper" ) .sublist ( "Eigensolver" );
    Teuchos::RCP<Teuchos::ParameterList> eigenListPtr =
        Teuchos::rcpFromRef ( ( eigenList ) );
    std::string eigenvaluesFileName =
        outputDirectory.string()  + "/" + outputList.get<string> ( "Eigenvalues file name" );
    std::string eigenstateFileNameAppendix =
        outputList.get<string> ( "Eigenstate file name appendix" );

    Teuchos::RCP<Ginla::IO::StatsWriter> eigenStatsWriter =
        Teuchos::rcp( new Ginla::IO::StatsWriter( eigenvaluesFileName ) );

    Teuchos::RCP<Ginla::IO::SaveEigenData> glEigenSaver =
        Teuchos::RCP<Ginla::IO::SaveEigenData> ( new Ginla::IO::SaveEigenData ( eigenListPtr,
                                                                                grid,
                                                                                glsystem->getKomplex(),
                                                                                eigenStatsWriter,
                                                                                stateWriter ) );

    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy =
        glEigenSaver;
    eigenList.set ( "Save Eigen Data Method", "User-Defined" );
    eigenList.set ( "User-Defined Save Eigen Data Name", "glSaveEigenDataStrategy" );
    eigenList.set ( "glSaveEigenDataStrategy", glSaveEigenDataStrategy );
#endif

    // ---------------------------------------------------------------------------
    // Create all possible Epetra_Operators.
    Teuchos::RCP<Epetra_RowMatrix> J = glsystem->getJacobian();
//  Teuchos::RCP<Epetra_RowMatrix> M = glsystem->getPreconditioner();

    // Create the linear system.
    // Use the TimeDependent interface for computation of shifted matrices.
    Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = glsystem;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glsystem;
//  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = glsystem;

    Teuchos::ParameterList& nlPrintParams =
        paramList->sublist ( "NOX" ) .sublist ( "Printing" );

    Teuchos::ParameterList& lsParams =
        paramList->sublist ( "NOX" ) .sublist ( "Direction" ) .sublist ( "Newton" ) .sublist ( "Linear Solver" );

//  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(
//      new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams, iJac, J, iPrec, M, *soln));

    NOX::Epetra::Vector cloneVector( Epetra_Vector( *glsystem->getExtendedMap() ) );
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
        Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( nlPrintParams,
                                                              lsParams,
                                                              iReq,
                                                              iJac,
                                                              J,
                                                              cloneVector ) );

    Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = glsystem;
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    
    // translate parameters into a LOCA list
    LOCA::ParameterVector locaParams =
          *(Ginla::Helpers::teuchosParameterList2locaParameterVector( glParameters ));

    // get initial guess
    double chi = 0.0;
    NOX::Epetra::Vector initialGuess ( glsystem->createInterfaceState( psi, chi ),
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
    // ---------------------------------------------------------------------------
    // Set up the NOX status tests
    Teuchos::ParameterList& noxList = paramList->sublist ( "NOX", true );
    double tol = noxList.get<double> ( "Tolerance" );
    int maxNonlinarSteps = noxList.get<int> ( "Max steps" );
    Teuchos::RCP<NOX::StatusTest::NormF> normF =
        Teuchos::rcp ( new NOX::StatusTest::NormF ( tol ) );
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxIters =
        Teuchos::rcp ( new NOX::StatusTest::MaxIters ( maxNonlinarSteps ) );
    Teuchos::RCP<NOX::StatusTest::Generic> comboOR =
        Teuchos::rcp ( new NOX::StatusTest::Combo ( NOX::StatusTest::Combo::OR, normF, maxIters ) );
    // ---------------------------------------------------------------------------
    // Set up the LOCA status tests
    Teuchos::RCP<LOCA::StatusTest::MaxIters> maxLocaStepsTest;
    try {
        maxLocaStepsTest = Teuchos::rcp ( new LOCA::StatusTest::MaxIters ( maxLocaSteps ) );
    }
    catch ( char const * e ) {
        std::cerr << e << std::endl;
        return 1;
    }
    // ---------------------------------------------------------------------------
    // Create the stepper
    Teuchos::RCP<LOCA::Thyra::Group> grp2 = Teuchos::null;
    Teuchos::RCP<LOCA::Stepper> stepper =
        Teuchos::rcp ( new LOCA::Stepper ( globalData,
                                           grp,
                                           maxLocaStepsTest,
                                           comboOR,
                                           paramList ) );
    // ---------------------------------------------------------------------------

    // make sure that the stepper starts off with the correct starting value

    // pass pointer to stepper to glsystem to be able to read stats from the stepper in there
    glsystem->setLocaStepper ( stepper );
#ifdef HAVE_LOCA_ANASAZI
    glEigenSaver->setLocaStepper ( stepper );
#endif

    // ---------------------------------------------------------------------------
    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status;
    try
    {
        status = stepper->run();
    }
    catch ( char const* e )
    {
        std::cerr << "Exception raised: " << e << std::endl;
    }
    // ---------------------------------------------------------------------------

    // clean up
    LOCA::destroyGlobalData ( globalData );
    glsystem->releaseLocaStepper();
#ifdef HAVE_LOCA_ANASAZI
    glEigenSaver->releaseLocaStepper();
#endif

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    // Final return value (0 = successful, non-zero = failure)
    return status;
}
// =============================================================================
