// see <http://old.nabble.com/Undefined-reference-to-%27main%27-with-Boost-Test.-Why--td15986217.html>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GRNN jacTest

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Vector.h>
#include <LOCA_Parameter_Vector.H>

#include "Ginla_StkMesh.hpp"
#include "Ginla_StkMeshReader.hpp"
#include "Ginla_MagneticVectorPotential.hpp"
#include "Ginla_JacobianOperator.hpp"

// ===========================================================================
// <http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/tutorials/hello-the-testing-world.html>
BOOST_AUTO_TEST_CASE( jac_test )
{
    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init ( &boost::unit_test::framework::master_test_suite().argc,
               &boost::unit_test::framework::master_test_suite().argv
             );
#endif

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm =
            Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm>  eComm =
            Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif
    // handle command line arguments
    Teuchos::CommandLineProcessor My_CLP;

    std::string inputFileName( "" );
    My_CLP.setOption ( "input", &inputFileName, "Input state file", true );

    // print warning for unrecognized arguments
    My_CLP.recogniseAllOptions ( true );

    // finally, parse the command line
    TEUCHOS_ASSERT_EQUALITY( My_CLP.parse ( boost::unit_test::framework::master_test_suite().argc,
                                            boost::unit_test::framework::master_test_suite().argv
                                          ),
                             Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL
                           );
    // =========================================================================
    // Read the data from the file.
    Teuchos::ParameterList data;
    Ginla::StkMeshRead( *eComm, inputFileName, data );

    // Cast the data into something more accessible.
    Teuchos::RCP<Ginla::StkMesh>     & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );
    Teuchos::RCP<Epetra_Vector>      & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );
    Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::ParameterList           & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

    Teuchos::RCP<Ginla::MagneticVectorPotential> mvp;
    double mu = 1.0e-2;
    mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential ( mesh, mvpValues, mu ) );

    Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
        Teuchos::rcp( new LOCA::ParameterVector() );
    mvpParameters->addParameter( "mu", mu );

    // create the jacobian operator
    Teuchos::RCP<Ginla::JacobianOperator> jac =
        Teuchos::rcp( new Ginla::JacobianOperator( mesh, thickness, mvp, z ) );

    // Create test vectors.
    // (a) vector os ones
    const Epetra_Map & map = jac->OperatorDomainMap();
    Teuchos::RCP<Epetra_Vector> s0 = Teuchos::rcp( new Epetra_Vector(map) );
    s0->PutScalar( 1.0 );
    Teuchos::RCP<Epetra_Vector> t0 = Teuchos::rcp( new Epetra_Vector(map) );
    jac->Apply( *s0, *t0 );

    // (b) [ 1, 0, 1, 0, ... ]
    Teuchos::RCP<Epetra_Vector> s1 = Teuchos::rcp( new Epetra_Vector(map) );
    Teuchos::RCP<Epetra_Vector> t1 = Teuchos::rcp( new Epetra_Vector(map) );
    double one  = 1.0;
    double zero = 0.0;
    for ( int k=0; k<map.NumMyPoints(); k++ )
    {
        if ( map.GID(k) % 2 )
            s1->ReplaceMyValues( 1, &one, &k );
        else
            s1->ReplaceMyValues( 1, &zero, &k );
    }
    jac->Apply( *s1, *t1 );

    // (b) [ 0, 1, 0, 1, ... ]
    Teuchos::RCP<Epetra_Vector> s2 = Teuchos::rcp( new Epetra_Vector(map) );
    Teuchos::RCP<Epetra_Vector> t2 = Teuchos::rcp( new Epetra_Vector(map) );
    for ( int k=0; k<map.NumMyPoints(); k++ )
    {
        if ( map.GID(k) % 2 )
            s2->ReplaceMyValues( 1, &zero, &k );
        else
            s2->ReplaceMyValues( 1, &one, &k );
    }
    jac->Apply( *s2, *t2 );

    // compute matrix norms as hashes
    double normT0;
    double normT1;
    double normT2;
    t0->Norm1( &normT0 );
    t1->Norm1( &normT1 );
    t2->Norm1( &normT2 );

    // control values
    double controlNormT0 = 20.2913968894543;
    double controlNormT1 = 0.289990630357596;
    double controlNormT2 = 20.2899906303576;

    // check the values
    BOOST_CHECK_SMALL( normT0-controlNormT0, 1.0e-12 );
    BOOST_CHECK_SMALL( normT1-controlNormT1, 1.0e-12 );
    BOOST_CHECK_SMALL( normT2-controlNormT2, 1.0e-12 );

    // clean up
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return;
}
// ============================================================================
