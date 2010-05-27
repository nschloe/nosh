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

#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_ModelEvaluatorInterface.H>

#include <NOX_Epetra_MatrixFree.H>
#include <NOX_Epetra_FiniteDifference.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Solver_Factory.H>

#include <NOX_StatusTest_Combo.H>
#include <NOX_StatusTest_NormF.H>
#include <NOX_StatusTest_NormUpdate.H>
#include <NOX_StatusTest_NormWRMS.H>
#include <NOX_StatusTest_MaxIters.H>
#include <NOX_StatusTest_FiniteValue.H>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <EpetraExt_RowMatrixOut.h>

#include <Piro_Epetra_LOCASolver.hpp>

#include <boost/filesystem.hpp>

// #include "Ginla_IO_SaveNewtonData.h"
#include "Ginla_IO_NoxObserver.h"
#include "Ginla_ModelEvaluator_Default.h"
#include "Recti_Grid_Reader.h"
#include "Ginla_MagneticVectorPotential_Centered.h"
#include "Ginla_Operator_BCCentral.h"
#include "Ginla_Operator_BCInner.h"
#include "Ginla_Operator_BCOuter.h"
#include "Ginla_IO_StateWriter.h"

#include <boost/filesystem.hpp>
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// ===========================================================================
// <http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/tutorials/hello-the-testing-world.htmlhh>
BOOST_AUTO_TEST_CASE( five_step_piro_test )
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

    // finally, parse the command line
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;

    parseReturn = My_CLP.parse ( boost::unit_test::framework::master_test_suite().argc,
                                boost::unit_test::framework::master_test_suite().argv );

    Teuchos::RCP<Teuchos::ParameterList> paramList =
        Teuchos::rcp ( new Teuchos::ParameterList );
    std::cout << "Reading parameter list from \"" << xmlInputFileName << "\"."
              << std::endl;

    Teuchos::updateParametersFromXmlFile ( xmlInputFileName, paramList.get() );
    // =========================================================================
    // extract data of the parameter list
    Teuchos::ParameterList& ioList = paramList->sublist ( "IO", true );

    boost::filesystem::path xmlPath =
        boost::filesystem::path ( xmlInputFileName ).branch_path();

    boost::filesystem::path inputGuessFile = ioList.get<string> ( "Input guess" );
    if ( !inputGuessFile.empty() && inputGuessFile.root_directory().empty() ) // if inputGuessFile is not empty and is a relative path
        inputGuessFile = xmlPath / inputGuessFile;

    const std::string outputFormat = ioList.get<string> ( "Output format" );
    
    // set default directory to be the directory of the XML file itself
    boost::filesystem::path outputDirectory = ioList.get<string> ( "Output directory" );
    if ( outputDirectory.root_directory().empty() )
        // outputDirectory is empty or is a relative directory.
        outputDirectory = xmlPath / outputDirectory;
    // =========================================================================

    Teuchos::RCP<Ginla::State>         initState;
    Teuchos::RCP<Recti::Grid::Uniform> grid;
    Teuchos::ParameterList             problemParameters;
    
    Recti::Grid::Reader::read ( Comm,
                                inputGuessFile.string(),
                                initState,
                                grid,
                                problemParameters );

    // possibly overwrite the parameters
    Teuchos::ParameterList & overwriteParamsList = paramList->sublist ( "Overwrite parameter list", true ); 
    bool overwriteParameters = overwriteParamsList.get<bool> ( "Overwrite parameters" );
    if ( overwriteParameters )
    {
        Teuchos::ParameterList & overwritePList = overwriteParamsList.sublist( "Parameters", true );
        problemParameters.setParameters( overwritePList );
        
        // possibly update the scaling of the grid
        grid->updateScaling( problemParameters.get<double>("scaling") );
    }
    
    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( problemParameters.get<double> ( "H0" ),
                                                                      problemParameters.get<double> ( "scaling" ) ) );

    Teuchos::RCP<Ginla::Komplex::LinearProblem> komplex =
        Teuchos::rcp( new Ginla::Komplex::LinearProblem( eComm, initState->getPsi()->getMap() ) );
        
    // setup the data output
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter =
        Teuchos::rcp( new Ginla::IO::StateWriter( outputDirectory.string(),
                                                  "solution",
                                                  "VTI",
                                                  1000 ) );
           
    // create the operator
    Teuchos::RCP<Ginla::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::Operator::BCCentral ( grid,
                                                        A,
                                                        initState->getPsi()->getMap(),
                                                        initState->getPsi()->getMap() ) );

    // Create the interface between NOX and the application
    // This object is derived from NOX::Epetra::Interface
    Teuchos::RCP<Ginla::ModelEvaluator::Default> glModel = 
              Teuchos::rcp(new Ginla::ModelEvaluator::Default( glOperator,
                                                               komplex,
                                                               *initState,
                                                               problemParameters ) );
                                                              
    Teuchos::RCP<Ginla::IO::NoxObserver> observer = Teuchos::null;

    Teuchos::RCP<Teuchos::ParameterList> piroParams =
        Teuchos::rcp(new Teuchos::ParameterList("Piro Parameters"));
    Teuchos::updateParametersFromXmlFile(xmlInputFileName, piroParams.get());

    // Use these two objects to construct a Piro solved application 
    //   EpetraExt::ModelEvaluator is  base class of all Piro::Epetra solvers
    Teuchos::RCP<EpetraExt::ModelEvaluator> piro
        = Teuchos::rcp( new Piro::Epetra::LOCASolver( piroParams,
                                                      glModel,
                                                      observer
                                                    )
                      );

    // Now the (somewhat cumbersome) setting of inputs and outputs
    EpetraExt::ModelEvaluator::InArgs inArgs = glModel->createInArgs();

    Teuchos::RCP<Epetra_Vector> p1 =
        Teuchos::rcp( new Epetra_Vector(*(piro->get_p_init(0))) );
//     int numParams = p1->MyLength(); // Number of parameters in p1 vector
    (*p1)[0] = problemParameters.get<double> ( "H0" );
    inArgs.set_p( 0, p1 );

    // Set output arguments to evalModel call
    EpetraExt::ModelEvaluator::OutArgs outArgs = piro->createOutArgs();

    // Solution vector is returned as extra response vector
    Teuchos::RCP<Epetra_Vector> gx = Teuchos::rcp(new Epetra_Vector(*(piro->get_g_map(0))));
    outArgs.set_g( 0, gx );
    
    // evaluate the model: go, go!
    piro->evalModel( inArgs, outArgs );
    
    // ------------------------------------------------------------------------
    // fetch the solution
    // outArgs.get_g(0) must be gx
    BOOST_ASSERT( !outArgs.get_g(0).is_null() ); 
    Teuchos::RCP<Ginla::State> solutionState = glModel->createState( *(outArgs.get_g(0)) );
    // ------------------------------------------------------------------------
    // read reference solution from file
    // For technical reasons, the reader can only accept ComplexMultiVectors.
    Teuchos::RCP<Ginla::State> refState;
    Recti::Grid::Reader::read ( Comm, expSolFileName, refState, grid, problemParameters );
    // ------------------------------------------------------------------------
    // compare the results:
    // get final solution
    
    solutionState->save( "test.vti" );
    
    Teuchos::RCP<Ginla::State> diff = solutionState;
    
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
