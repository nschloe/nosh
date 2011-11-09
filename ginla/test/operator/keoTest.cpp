// see <http://old.nabble.com/Undefined-reference-to-%27main%27-with-Boost-Test.-Why--td15986217.html>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GRNN keoTest

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
#include "Ginla_KeoFactory.hpp"

// ===========================================================================
// <http://www.boost.org/doc/libs/1_42_0/libs/test/doc/html/tutorials/hello-the-testing-world.html>
BOOST_AUTO_TEST_CASE( keo_test )
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

    Teuchos::RCP<Ginla::KeoFactory> keoFactory =
        Teuchos::rcp( new Ginla::KeoFactory( mesh, thickness, mvp ) );

    Teuchos::RCP<Epetra_FECrsGraph> keoGraph =
        Teuchos::rcp( new Epetra_FECrsGraph( keoFactory->buildKeoGraph() ) );

    // create the kinetic energy operator
    Teuchos::RCP<Epetra_FECrsMatrix> keoMatrix;
    keoMatrix = Teuchos::rcp( new Epetra_FECrsMatrix( Copy, *keoGraph ) );
    keoFactory->updateParameters( mvpParameters );
    keoFactory->buildKeo( *keoMatrix );

    // Make sure the matrix is indeed positive definite, and not
    // negative definite. Belos needs that (2010-11-05).
    keoMatrix->Scale( -1.0 );

    // Compute matrix norms as hashes.
    // Don't check for NormFrobenius() as this one doesn't work for matrices
    // with overlapping maps.
    double normOne = keoMatrix->NormOne();
    double normInf = keoMatrix->NormInf();

    // control values
    double controlNormOne = 10.1456474156918;
    double controlNormInf = 10.1456474156918;
    
    // check the values
    BOOST_CHECK_SMALL( normOne-controlNormOne, 1.0e-13 );
    BOOST_CHECK_SMALL( normInf-controlNormInf, 1.0e-13 );

    const Epetra_Map & map = keoMatrix->DomainMap();

    // Add up all the entries of the matrix.
    Epetra_Vector e( map );
    e.PutScalar( 1.0 );
    Epetra_Vector Ke( map );
    keoMatrix->Apply( e, Ke );
    double sum[1];
    e.Dot( Ke, sum );
    double controlSum = 0.00844428504187249;
    BOOST_CHECK_SMALL( sum[0]-controlSum, 1.0e-13 );

    // Sum over all the "real parts" of the matrix.
    // Remember that a 2x2 block corresponding to z is composed as
    // [ Re(z) -Im(z) ]
    // [ Im(z)  Re(z) ].
    // Build vector [ 1, 0, 1, 0, ... ]:
    double one  = 1.0;
    double zero = 0.0;
    Epetra_Vector s0( map );
    for ( int k=0; k<map.NumMyPoints(); k++ )
    {
        if ( map.GID(k) % 2 )
            s0.ReplaceMyValues( 1, &one, &k );
        else
            s0.ReplaceMyValues( 1, &zero, &k );
    }
    Epetra_Vector t0( map );
    keoMatrix->Apply( s0, t0 );
    s0.Dot( t0, sum );
    double controlSumReal = 0.0042221425209367988;
    BOOST_CHECK_SMALL( sum[0]-controlSumReal, 1.0e-13 );

    // Sum over all the "imaginary parts" of the matrix.
    Epetra_Vector s1( map );
    for ( int k=0; k<map.NumMyPoints(); k++ )
    {
        if ( map.GID(k) % 2 )
            s1.ReplaceMyValues( 1, &zero, &k );
        else
            s1.ReplaceMyValues( 1, &one, &k );
    }
    Epetra_Vector t1( map );
    keoMatrix->Apply( s0, t1 );
    s1.Dot( t0, sum );
    double controlSumImag = 0.0; // Matrix is Hermitian
    BOOST_CHECK_SMALL( sum[0]-controlSumImag, 1.0e-13 );

    // clean up
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return;
}
// ============================================================================
