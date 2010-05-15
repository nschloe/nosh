// see <http://old.nabble.com/Undefined-reference-to-%27main%27-with-Boost-Test.-Why--td15986217.html>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GRNN vorticityTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <LOCA_Epetra_Factory.H>
#include <LOCA_Epetra_Group.H>
#include <LOCA_StatusTest_MaxIters.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <NOX_StatusTest_NormF.H>
#include <NOX_StatusTest_MaxIters.H>
#include <NOX_StatusTest_Combo.H>
// #include <NOX_Epetra_Group.H>

#include "Ginla_State.h"

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
BOOST_AUTO_TEST_CASE( vorticity_test )
{
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( NULL, NULL );
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
    boost::filesystem::path xmlPath =
            boost::filesystem::path ( xmlInputFileName ).branch_path();

    Teuchos::RCP<Teuchos::ParameterList> paramList =
            Teuchos::rcp ( new Teuchos::ParameterList );            
    Teuchos::updateParametersFromXmlFile ( xmlInputFileName,
                                           paramList.get()
                                         );
    // ------------------------------------------------------------------------
    // iterate through the list of example solutions, calculate their energies,
    // and compare with the given value
    Teuchos::RCP<ComplexMultiVector> psiM;
    Teuchos::RCP<ComplexVector> psi;
    Teuchos::RCP<Recti::Grid::Uniform> grid;
    Teuchos::ParameterList glParameters;
    for ( Teuchos::ParameterList::ConstIterator k=paramList->begin(); k!=paramList->end(); ++k )
    {
        std::string label = paramList->name( k );   
        Teuchos::ParameterList & stateSublist = paramList->sublist( label );

        boost::filesystem::path filePath = stateSublist.get<string> ( "File name" );
        filePath = xmlPath / filePath;
        
        int vorticityRef = stateSublist.get<int>( "Vorticity" );
        
        // read the state
        Recti::Grid::Reader::read ( Comm, filePath.string(), psiM, grid, glParameters );
        TEUCHOS_ASSERT_EQUALITY ( psiM->getNumVectors(), 1 );
        psi = psiM->getVectorNonConst ( 0 );
        
        Teuchos::RCP<Ginla::State> state = 
            Teuchos::rcp( new Ginla::State( psiM->getVectorNonConst(0),
                                            grid ) );       
        double vorticity = state->getVorticity();

        BOOST_CHECK_EQUAL( vorticity, vorticityRef );
    }
    // ------------------------------------------------------------------------
    // clean up
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    // ------------------------------------------------------------------------

    return;
}
// ============================================================================
