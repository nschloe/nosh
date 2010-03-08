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

#include "ginzburgLandau.h"
#include "GL_Operator_BCCentral.h"
#include "GL_Operator_BCInner.h"

#include "Recti_Domain_Square.h"

#include "GL_LocaSystem_Bordered.h"
#include "GL_IO_SaveEigenData.h"

#include "GL_Helpers.h"
#include "GL_StatsWriter.h"

#include "Recti_Grid_Reader.h"


// =============================================================================
unsigned int
numDigits ( const int i )
{
    int numDigits = 0;
    int ii = i;
    if ( ii < 0 )
        ii = -ii;

    while ( ii > 0 )
    {
        numDigits++;
        ii/=10;
    }
    return numDigits;
}
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
    std::string contFileFormat =
        outputList.get<string> ( "Continuation file format" );
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
    boost::filesystem::path inputGuessFile = initialGuessList.get<string> ( "Solution" );
    if ( !inputGuessFile.empty() && inputGuessFile.root_directory().empty() ) // if inputGuessFile is a relative path
        inputGuessFile = xmlPath / inputGuessFile;

    if ( !inputGuessFile.empty() )
    {
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

        if ( initialGuessList.isParameter ( "Nx" ) )
        {
            std::cerr << "Warning: Parameter 'Nx' *and* input guess file given."
                      << "Using value as in input guess file, discarding 'Nx'."
                      << std::endl;
        }
        if ( initialGuessList.isParameter ( "scaling" ) )
        {
            std::cerr << "Warning: Parameter 'scaling' *and* input guess file given."
                      << "Using value as in input guess file, discarding 'scaling'."
                      << std::endl;
        }
        if ( initialGuessList.isParameter ( "H0" ) )
        {
            std::cerr << "Warning: Parameter 'H0' *and* input guess file given."
                      << "Using value as in input guess file, discarding 'H0'."
                      << std::endl;
        }

    }
    else // no or empty input guess file
    {
        double scaling = initialGuessList.get<double> ( "scaling" );
        double H0      = initialGuessList.get<double> ( "H0" );

        glParameters.set ( "scaling", scaling );
        glParameters.set ( "H0", H0 );

        double edgeLength = 1.0;
        Teuchos::RCP<Recti::Domain::Abstract> domain = Teuchos::rcp ( new Recti::Domain::Square ( edgeLength ) );

        int    Nx = initialGuessList.get<int> ( "Nx" );
        double h  = edgeLength / Nx;
        Teuchos::RCP<Recti::Grid::Uniform> grid = Teuchos::rcp ( new Recti::Grid::Uniform ( domain, h ) );
        grid->updateScaling ( scaling );

        // set initial guess
        int numComplexUnknowns = grid->getNumGridPoints();
        Teuchos::RCP<Tpetra::Map<Thyra::Ordinal> > complexMap =
            Teuchos::rcp ( new Tpetra::Map<Thyra::Ordinal> ( numComplexUnknowns, 0, Comm ) );

        psi = Teuchos::rcp ( new ComplexVector ( complexMap ) );
        double_complex alpha ( 1.0, 0.0 );
        psi->putScalar ( alpha );
    }
    // ---------------------------------------------------------------------------

    double h0      = glParameters.get<double> ( "H0" );
    double scaling = glParameters.get<double> ( "scaling" );
    Teuchos::RCP<GL::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new GL::MagneticVectorPotential::Centered ( h0, scaling ) );

    // create the operator
    Teuchos::RCP<GL::Operator::Virtual> glOperator =
        Teuchos::rcp ( new GL::Operator::BCCentral ( grid, A ) );

        
    std::string fn = outputDirectory.string() + "/" + contDataFileName;
    Teuchos::RCP<GL::StatsWriter> statsWriter =
        Teuchos::rcp( new GL::StatsWriter( fn ) );
    GinzburgLandau glProblem = GinzburgLandau ( glOperator,
                                                statsWriter,
                                                contFileFormat );

    Teuchos::RCP<GL::LocaSystem::Bordered> glsystem;

    Teuchos::ParameterList & stepperList = paramList->sublist ( "LOCA" ).sublist ( "Stepper" );
    int maxLocaSteps = stepperList.get<int> ( "Max Steps" );
    try
    {
        glsystem = Teuchos::rcp ( new GL::LocaSystem::Bordered ( glProblem,
                                  eComm,
                                  psi,
                                  outputDirectory.string(),
                                  contDataFileName,
                                  contFileBaseName,
                                  "",
                                  numDigits ( maxLocaSteps ) ) );
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

    // get the initial solution
    Teuchos::RCP<Epetra_Vector> soln = glsystem->getSolution();
    // ---------------------------------------------------------------------------

#ifdef HAVE_LOCA_ANASAZI
    Teuchos::ParameterList& eigenList =
        paramList->sublist ( "LOCA" ).sublist ( "Stepper" ) .sublist ( "Eigensolver" );
    Teuchos::RCP<Teuchos::ParameterList> eigenListPtr =
        Teuchos::rcpFromRef ( ( eigenList ) );
    std::string eigenvaluesFileName =
        outputList.get<string> ( "Eigenvalues file name" );
    std::string eigenstateFileNameAppendix =
        outputList.get<string> ( "Eigenstate file name appendix" );
    Teuchos::RCP<GL::IO::SaveEigenData> glEigenSaver =
        Teuchos::RCP<GL::IO::SaveEigenData> ( new GL::IO::SaveEigenData ( eigenListPtr, outputDirectory.string(),
                                   eigenvaluesFileName, contFileBaseName,
                                   eigenstateFileNameAppendix, glsystem,
                                   numDigits ( maxLocaSteps ) ) );

    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy =
        glEigenSaver;
    eigenList.set ( "glSaveEigenDataStrategy", glSaveEigenDataStrategy );
#endif

    // ---------------------------------------------------------------------------
    // Create all possible Epetra_Operators.
    Teuchos::RCP<Epetra_RowMatrix> J = glsystem->getJacobian();

    // Create the linear system.
    // Use the TimeDependent interface for computation of shifted matrices.
    Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = glsystem;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = glsystem;

    Teuchos::ParameterList& nlPrintParams =
        paramList->sublist ( "NOX" ) .sublist ( "Printing" );

    Teuchos::ParameterList& lsParams =
        paramList->sublist ( "NOX" ) .sublist ( "Direction" ) .sublist ( "Newton" ) .sublist ( "Linear Solver" );

    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
        Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( nlPrintParams, lsParams, iReq, iJac, J, *soln ) );

    Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = glsystem;
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.

    // translate parameters into a LOCA list
    LOCA::ParameterVector locaParams =
          *(GL::Helpers::teuchosParameterList2locaParameterVector( glParameters ));

    NOX::Epetra::Vector initialGuess ( soln, NOX::Epetra::Vector::CreateView );

    Teuchos::RCP<LOCA::Epetra::Group> grp =
        Teuchos::rcp ( new LOCA::Epetra::Group ( globalData, nlPrintParams, iTime,
                       initialGuess, linSys, linSys, locaParams ) );

    grp->setParams ( locaParams );
    // ---------------------------------------------------------------------------


    // ---------------------------------------------------------------------------
    // Get the vector from the Problem
    Teuchos::RCP<NOX::Epetra::Vector> noxSoln =
        Teuchos::rcp ( new NOX::Epetra::Vector ( soln, NOX::Epetra::Vector::CreateView ) );
    // ---------------------------------------------------------------------------


    // ---------------------------------------------------------------------------
    // Set up the NOX status tests
    // ---------------------------------------------------------------------------
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

    // ---------------------------------------------------------------------------
    // Set up the LOCA status tests
    // ---------------------------------------------------------------------------
    Teuchos::RCP<LOCA::StatusTest::MaxIters> maxLocaStepsTest;
    try {
        maxLocaStepsTest = Teuchos::rcp ( new LOCA::StatusTest::MaxIters ( maxLocaSteps ) );
    }
    catch ( char const * e ) {
        std::cerr << e << std::endl;
        return 1;
    }
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // read in initial null vector and convert it to a glsystem-compliant vector
    boost::filesystem::path initialNullVectorFile = initialGuessList.get<string> ( "Null vector" );
    if ( !initialNullVectorFile.empty() && initialNullVectorFile.root_directory().empty() ) // if initialNullVectorFile is a relative path
        initialNullVectorFile = xmlPath / initialNullVectorFile;
    
    Teuchos::RCP<ComplexVector> initialNullVector;
    try
    {
        Teuchos::RCP<ComplexMultiVector> initialNullVectorM;
        Recti::Grid::Reader::read ( Comm, initialNullVectorFile.string(), initialNullVectorM, grid, glParameters );
        TEUCHOS_ASSERT_EQUALITY ( initialNullVectorM->getNumVectors(), 1 );
        initialNullVector = initialNullVectorM->getVectorNonConst ( 0 );
    }
    catch ( std::exception & e )
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch ( ... )
    {
        std::cerr << "Unknown exception caught." << std::endl;
        return 1;
    }
    // convert the complex vector into a (real-valued) GlSystem-compliant vector
    Teuchos::RCP<Epetra_Vector> glsystemInitialNullVector =
        glsystem->getGlKomplex()->complex2real ( initialNullVector );
    // ---------------------------------------------------------------------------


    // ---------------------------------------------------------------------------
    // add LOCA options which cannot be provided in the XML file
    Teuchos::ParameterList & bifList =
        paramList->sublist ( "LOCA" ).sublist ( "Bifurcation" );
        
    // create extended vectors
    Teuchos::RCP<const Epetra_BlockMap> regularMap = Teuchos::rcpFromRef ( glsystemInitialNullVector->Map() );
    Teuchos::RCP<const Epetra_Map> extendedMap = glsystem->getExtendedMap();
    Teuchos::RCP<Epetra_Vector> glsystemInitialNullVectorExt =
        Teuchos::rcp( new Epetra_Vector( *extendedMap, true ) );
    Teuchos::RCP<Epetra_Import> importFromRegularMap =
         Teuchos::rcp (new Epetra_Import(*extendedMap,*regularMap)); 
    glsystemInitialNullVectorExt->Import( *glsystemInitialNullVector,*importFromRegularMap,Insert );
    
        
    Teuchos::RCP<NOX::Abstract::Vector> lengthNormVec =
        Teuchos::rcp ( new NOX::Epetra::Vector ( *glsystemInitialNullVectorExt ) );
//     lengthNormVec->init(1.0);
    
    bifList.set ( "Length Normalization Vector", lengthNormVec );

    Teuchos::RCP<NOX::Abstract::Vector> initialNullAbstractVec =
        Teuchos::rcp ( new NOX::Epetra::Vector ( *glsystemInitialNullVectorExt ) );
//     initialNullAbstractVec->init(1.0);
    bifList.set ( "Initial Null Vector", initialNullAbstractVec );
    // ---------------------------------------------------------------------------


    // ---------------------------------------------------------------------------
    // Create the stepper
    Teuchos::RCP<LOCA::Stepper> stepper;
    try {
        stepper = Teuchos::rcp ( new LOCA::Stepper ( globalData,
                                                     grp,
                                                     maxLocaStepsTest,
                                                     comboOR,
                                                     paramList ) );
    }
    catch (...) {
        std::cerr << "Exception caught. Abort." << std::endl;
        return 1;
    }
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
