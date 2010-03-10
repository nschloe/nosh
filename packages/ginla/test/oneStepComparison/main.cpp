// see <http://old.nabble.com/Undefined-reference-to-%27main%27-with-Boost-Test.-Why--td15986217.html>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GRNN MyTest

#include <boost/test/unit_test.hpp>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <boost/filesystem.hpp>
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

BOOST_AUTO_TEST_SUITE( my_suite )
// ===========================================================================
// <http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/tutorials/hello-the-testing-world.htmlhh>
BOOST_AUTO_TEST_CASE( my_test )
{
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( &argc,&argv );
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
    // Read the XML file
    std::string xmlInputFileName = "data/conf.xml";
    Teuchos::RCP<Teuchos::ParameterList> paramList =
            Teuchos::rcp ( new Teuchos::ParameterList );
    try
    {
        Teuchos::updateParametersFromXmlFile ( xmlInputFileName, paramList.get() );
    }
    catch ( std::exception &e )
    {
        std::cerr << e.what() << "\n"
                  << "Check your XML file for syntax errors." << std::endl;
        return EXIT_FAILURE;
    }
    // ------------------------------------------------------------------------
    // extract data for the output the parameter list
    Teuchos::ParameterList outputList;
    try
    {
        outputList = paramList->sublist ( "Output", true );
    }
    catch ( std::exception &e )
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
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
    // ------------------------------------------------------------------------
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
    // ------------------------------------------------------------------------

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
    Teuchos::RCP<Ginla::StatsWriter> statsWriter = 
        Teuchos::rcp( new Ginla::StatsWriter( fn ) );

    GinzburgLandau glProblem = GinzburgLandau ( glOperator,
                                                statsWriter,
                                                outputFormat
                                              );

    Teuchos::RCP<Ginla::LocaSystem::Bordered> glsystem;

    Teuchos::ParameterList & stepperList = paramList->sublist ( "LOCA" ).sublist ( "Stepper" );
    int maxLocaSteps = stepperList.get<int> ( "Max Steps" );

    try
    {
        glsystem = Teuchos::rcp ( new Ginla::LocaSystem::Bordered ( glProblem,
                                  eComm,
                                  psi,
                                  outputDirectory.string(),
                                  contDataFileName,
                                  contFileBaseName,
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

    // ------------------------------------------------------------------------
    // Create the necessary objects
    // ------------------------------------------------------------------------
    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
        Teuchos::rcp ( new LOCA::Epetra::Factory );

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
        LOCA::createGlobalData ( paramList, epetraFactory );

    // get the initial solution
    Teuchos::RCP<Epetra_Vector> soln = glsystem->getSolution();
    // ------------------------------------------------------------------------

#ifdef HAVE_LOCA_ANASAZI
    Teuchos::ParameterList& eigenList =
        paramList->sublist ( "LOCA" ).sublist ( "Stepper" ) .sublist ( "Eigensolver" );
    Teuchos::RCP<Teuchos::ParameterList> eigenListPtr =
        Teuchos::rcpFromRef ( ( eigenList ) );
    std::string eigenvaluesFileName =
        outputList.get<string> ( "Eigenvalues file name" );
    std::string eigenstateFileNameAppendix =
        outputList.get<string> ( "Eigenstate file name appendix" );
    Teuchos::RCP<Ginla::IO::SaveEigenData> glEigenSaver =
        Teuchos::RCP<Ginla::IO::SaveEigenData> ( new Ginla::IO::SaveEigenData ( eigenListPtr, outputDirectory.string(),
                                   eigenvaluesFileName, contFileBaseName,
                                   eigenstateFileNameAppendix, glsystem,
                                   numDigits ( maxLocaSteps ) ) );

    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy> glSaveEigenDataStrategy =
        glEigenSaver;
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
        paramList->sublist ( "NOX" ) .sublist ( "Printing" );

    Teuchos::ParameterList& lsParams =
        paramList->sublist ( "NOX" ) .sublist ( "Direction" ) .sublist ( "Newton" ) .sublist ( "Linear Solver" );

//  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = Teuchos::rcp(
//      new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams, iJac, J, iPrec, M, *soln));


    std::cout << J.is_null() << std::endl;
    if ( !soln.is_valid_ptr() || soln.is_null() )
    {
      std::cout << "soln not properly initialized. Abort." << std::endl;
      return 1;
    }
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
        Teuchos::rcp ( new NOX::Epetra::LinearSystemAztecOO ( nlPrintParams, lsParams, iReq, iJac, J, *soln ) );

    Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = glsystem;
    // ------------------------------------------------------------------------
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    
    // translate parameters into a LOCA list
    LOCA::ParameterVector locaParams =
          *(Ginla::Helpers::teuchosParameterList2locaParameterVector( glParameters ));

    NOX::Epetra::Vector initialGuess ( soln, NOX::Epetra::Vector::CreateView );

    Teuchos::RCP<LOCA::Epetra::Group> grp =
        Teuchos::rcp ( new LOCA::Epetra::Group ( globalData, nlPrintParams, iTime,
                       initialGuess, linSys, linSys, locaParams ) );

    grp->setParams ( locaParams );
    // ------------------------------------------------------------------------
    // Get the vector from the Problem
    Teuchos::RCP<NOX::Epetra::Vector> noxSoln =
        Teuchos::rcp ( new NOX::Epetra::Vector ( soln, NOX::Epetra::Vector::CreateView ) );
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
        Teuchos::rcp ( new NOX::StatusTest::Combo ( NOX::StatusTest::Combo::OR, normF, maxIters ) );
    // ------------------------------------------------------------------------
    // Set up the LOCA status tests
    Teuchos::RCP<LOCA::StatusTest::MaxIters> maxLocaStepsTest;
    try {
        maxLocaStepsTest = Teuchos::rcp ( new LOCA::StatusTest::MaxIters ( maxLocaSteps ) );
    }
    catch ( char const * e ) {
        std::cerr << e << std::endl;
        return 1;
    }
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
    try
    {
        status = stepper->run();
    }
    catch ( char const* e )
    {
        std::cerr << "Exception raised: " << e << std::endl;
    }
    // ------------------------------------------------------------------------
    // retrieve solution
    const NOX::Epetra::Group & finalGroup =
        dynamic_cast<const NOX::Epetra::Group&> ( stepper->getSolutionGroup() );

    const Epetra_Vector & finalSolution =
        ( dynamic_cast<const NOX::Epetra::Vector&> ( finalGroup.getX() ) ).getEpetraVector();
    std::cout << "Fff" << std::endl;
    // ------------------------------------------------------------------------
    // read expected solution from file
    
    // ------------------------------------------------------------------------
    // compare the results
    bool status = false;
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
    BOOST_CHECK( status );
}
// ============================================================================
BOOST_AUTO_TEST_CASE( my_other_test )
{
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( &argc,&argv );
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

    bool status = true;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    BOOST_CHECK( status );
}
// ============================================================================
BOOST_AUTO_TEST_SUITE_END()