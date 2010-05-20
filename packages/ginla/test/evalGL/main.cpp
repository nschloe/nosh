// see <http://old.nabble.com/Undefined-reference-to-%27main%27-with-Boost-Test.-Why--td15986217.html>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GRNN evalGl

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include "Ginla_MagneticVectorPotential_Centered.h"
#include "Ginla_Operator_BCCentral.h"
#include "Ginla_LocaSystem_Default.h"

#include "Recti_Grid_Reader.h"

#include <Epetra_CrsMatrix.h>
#include <EpetraExt_HDF5.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <boost/filesystem.hpp>
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

BOOST_AUTO_TEST_SUITE( eval_gl_test )
// ===========================================================================
// <http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/tutorials/hello-the-testing-world.htmlhh>
BOOST_AUTO_TEST_CASE( eval_operator_test )
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

    std::string inputStateFile = "";
    My_CLP.setOption ( "input-state-file", &inputStateFile,
                       "VTK/VTI file containing the input state", true );
    
    std::string expSolFileName = "";
    My_CLP.setOption ( "reference-residual-file", &expSolFileName,
                       "VTK/VTI file containing the expected residual", true );

    std::string expJacobianFileName = "";
    My_CLP.setOption ( "reference-jacobian-file", &expJacobianFileName,
                       "HDF5 file containing the expected Jacobian", true );
                       
    // print warning for unrecognized arguments
    My_CLP.recogniseAllOptions ( true );

    // don't throw exceptions
    My_CLP.throwExceptions ( true );

    // finally, parse the stuff!
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
    
    parseReturn = My_CLP.parse ( boost::unit_test::framework::master_test_suite().argc,
                                 boost::unit_test::framework::master_test_suite().argv );
    // ------------------------------------------------------------------------
    // read the input state
    Teuchos::RCP<Ginla::State>         state;
    Teuchos::RCP<Recti::Grid::Uniform> grid;
    Teuchos::ParameterList             glParameters;
    Recti::Grid::Reader::read ( Comm, inputStateFile, state, grid, glParameters );
    // ------------------------------------------------------------------------
    // set the operator
    double h0 = 0.9;
    double scaling = 7.0;
    
    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( h0, scaling ) );

    // create the operator
    Teuchos::RCP<const ComplexMap> map = state->getPsi()->getMap();
    Teuchos::RCP<Ginla::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::Operator::BCCentral ( grid, A, map, map ) );
    // ------------------------------------------------------------------------
    // evaluate
    Teuchos::RCP<Ginla::State> residual = glOperator->getF( state );
    // ------------------------------------------------------------------------    
    // read the reference state
    Teuchos::RCP<Ginla::State> referenceResidual;
    Recti::Grid::Reader::read ( Comm, expSolFileName, referenceResidual, grid, glParameters );
    // ------------------------------------------------------------------------
    // and compare
    residual->update( -1.0, *referenceResidual, 1.0 );
    
    double nrm = residual->normalizedScaledL2Norm();
    
    BOOST_CHECK_SMALL( nrm, 1.0e-10 );  
    // ------------------------------------------------------------------------
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return;
}
// ============================================================================
BOOST_AUTO_TEST_CASE( eval_locasystem_test )
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

    std::string inputStateFile = "";
    My_CLP.setOption ( "input-state-file", &inputStateFile,
                       "VTK/VTI file containing the input state", true );
    
    std::string expSolFileName = "";
    My_CLP.setOption ( "reference-residual-file", &expSolFileName,
                       "VTK/VTI file containing the expected residual", true );
                       
    std::string referenceJacobianFileName = "";
    My_CLP.setOption ( "reference-jacobian-file", &referenceJacobianFileName,
                       "HDF5 file containing the expected Jacobian", true );

    // print warning for unrecognized arguments
    My_CLP.recogniseAllOptions ( true );

    // don't throw exceptions
    My_CLP.throwExceptions ( true );

    // finally, parse the stuff!
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
    
    parseReturn = My_CLP.parse ( boost::unit_test::framework::master_test_suite().argc,
                                 boost::unit_test::framework::master_test_suite().argv );
    // ------------------------------------------------------------------------
    // read the input state
    Teuchos::RCP<Ginla::State>         state;
    Teuchos::RCP<Recti::Grid::Uniform> grid;
    Teuchos::ParameterList             glParameters;
    Recti::Grid::Reader::read ( Comm, inputStateFile, state, grid, glParameters );
    // ------------------------------------------------------------------------
    // set the operator
    double h0 = 0.9;
    double scaling = 7.0;
    
    Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A =
        Teuchos::rcp ( new Ginla::MagneticVectorPotential::Centered ( h0, scaling ) );

    // create the operator
    Teuchos::RCP<const ComplexMap> map = state->getPsi()->getMap();
    Teuchos::RCP<Ginla::Operator::Virtual> glOperator =
        Teuchos::rcp ( new Ginla::Operator::BCCentral ( grid, A, map, map ) );
    // ------------------------------------------------------------------------
    // create some dummy writers
    Teuchos::RCP<Ginla::IO::StatsWriter> statsWriter;
    Teuchos::RCP<Ginla::IO::StateWriter> stateWriter;
    // ------------------------------------------------------------------------
    // create the system
    Teuchos::RCP<Ginla::LocaSystem::Default> system
        = Teuchos::rcp( new Ginla::LocaSystem::Default ( glOperator,
                                                         eComm,
                                                         map,
                                                         statsWriter,
                                                         stateWriter ) );
    // ------------------------------------------------------------------------
    // evaluate the system in 'state'
    Teuchos::RCP<Epetra_Vector> x = system->createSystemVector( *state );
    Teuchos::RCP<Epetra_Vector> F = system->createSystemVector( *state );

    system->computeF( *x, *F );
    
    Teuchos::RCP<Ginla::State> residual = system->createState( *F );
    // ------------------------------------------------------------------------    
    // read the reference state
    Teuchos::RCP<Ginla::State> referenceResidual;
    Recti::Grid::Reader::read ( Comm, expSolFileName, referenceResidual, grid, glParameters );
    // ------------------------------------------------------------------------
    // and compare
    residual->update( -1.0, *referenceResidual, 1.0 );
    
    double nrm = residual->normalizedScaledL2Norm();
    
    BOOST_CHECK_SMALL( nrm, 1.0e-10 );  
    // ------------------------------------------------------------------------
    Teuchos::RCP<Epetra_CrsMatrix> jacobian = system->getJacobian();
    system->computeJacobian( *x, *jacobian );
    jacobian = system->getJacobian();
    // ------------------------------------------------------------------------
    // read the reference jacobian
    EpetraExt::HDF5 hdf5( *eComm );
    hdf5.Open( referenceJacobianFileName );
    Epetra_CrsMatrix * referenceJacobian = 0;
    std::string name = "jacobian";
    hdf5.Read( name, referenceJacobian );
    hdf5.Close();
    // ------------------------------------------------------------------------
    // compare the jacobians
    int numRowEntries, refNumRowEntries;
    double * values;
    double * refValues;
    int * indices;
    int * refIndices;
    for ( int row=0; row<referenceJacobian->NumMyRows(); row++ )
    {
       // extract the row information
       jacobian->ExtractMyRowView( row, numRowEntries, values, indices );
       referenceJacobian->ExtractMyRowView( row, refNumRowEntries, refValues, refIndices );
       // compare the entries
       BOOST_REQUIRE_EQUAL( numRowEntries, refNumRowEntries );
       for ( int k=0; k<numRowEntries; k++ )
       {
         BOOST_REQUIRE_EQUAL( indices[k], refIndices[k] );
         BOOST_REQUIRE_SMALL( values[k] - refValues[k], 1.0e-10 );
       }
    }
    // ------------------------------------------------------------------------    
//     // store the jacobian in a file
//     EpetraExt::HDF5 hdf5( *eComm );
//     hdf5.Create( inputStateFile + ".h5" );
//     hdf5.Write( "jacobian", *jac );
//     hdf5.Close();
    // ------------------------------------------------------------------------
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return;
}
// ============================================================================
BOOST_AUTO_TEST_SUITE_END()