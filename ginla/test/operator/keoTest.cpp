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
    Teuchos::ParameterList              problemParameters;
    Teuchos::RCP<Epetra_Vector>         z = Teuchos::null;
    Teuchos::RCP<Epetra_MultiVector>    mvpValues = Teuchos::null;
    Teuchos::RCP<Epetra_Vector>         thickness = Teuchos::null;
    Teuchos::RCP<Ginla::StkMesh> mesh = Teuchos::null;

    // read the state
    Ginla::StkMeshRead( *eComm,
                        inputFileName,
                        z,
                        mvpValues,
                        thickness,
                        mesh,
                        problemParameters
                    );

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

    // compute matrix norms as hashes
    double normOne = keoMatrix->NormOne();
    double normInf = keoMatrix->NormInf();
    double normFro = keoMatrix->NormFrobenius();

    // control values
    double controlNormOne = 60.0;
    double controlNormInf = 60.0;
    double controlNormFro = 71.414284285428437;
    
    // check the values
    BOOST_CHECK_SMALL( normOne-controlNormOne, 1.0e-13 );
    BOOST_CHECK_SMALL( normInf-controlNormInf, 1.0e-13 );
    BOOST_CHECK_SMALL( normFro-controlNormFro, 1.0e-13 );

    // clean up
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return;
}
// ============================================================================
