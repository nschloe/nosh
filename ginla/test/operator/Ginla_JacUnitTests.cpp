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
#include "Ginla_JacobianOperator.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace {

// =============================================================================
void
testJac( const std::string & inputFileNameBase,
         const double mu,
         const double controlNormT0,
         const double controlNormT1,
         const double controlNormT2,
         Teuchos::FancyOStream & out,
         bool & success )
{
    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Teuchos::RCP<Epetra_MpiComm> eComm =
            Teuchos::rcp<Epetra_MpiComm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    Teuchos::RCP<Epetra_SerialComm> eComm =
            Teuchos::rcp<Epetra_SerialComm> ( new Epetra_SerialComm() );
#endif

    std::string inputFileName;
    if ( eComm->NumProc() == 1 )
        inputFileName = inputFileNameBase + ".e";
    else
        inputFileName = inputFileNameBase + "-balanced.par";
    // =========================================================================
    // Read the data from the file.
    Teuchos::ParameterList data;
    Ginla::StkMeshRead( *eComm, inputFileName, data );

    // Cast the data into something more accessible.
    Teuchos::RCP<Ginla::StkMesh>     & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );
    Teuchos::RCP<Epetra_Vector>      & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::RCP<const Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<const Epetra_MultiVector>() );
    Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::ParameterList           & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

    Teuchos::RCP<Ginla::MagneticVectorPotential> mvp;
    mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential ( mesh, mvpValues, mu ) );

    Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
        Teuchos::rcp( new LOCA::ParameterVector() );
    mvpParameters->addParameter( "mu", mu );

    // create a keo factory
    Teuchos::RCP<Ginla::KeoFactory> keoFactory =
        Teuchos::rcp( new Ginla::KeoFactory( mesh, thickness, mvp ) );

    // create the jacobian operator
    Teuchos::RCP<Ginla::JacobianOperator> jac =
        Teuchos::rcp( new Ginla::JacobianOperator( mesh, thickness, keoFactory, z ) );

    // Create test vectors.
    // (a) vector os ones
    const Epetra_Map & map = jac->OperatorDomainMap();
    Teuchos::RCP<Epetra_Vector> s0 = Teuchos::rcp( new Epetra_Vector(map) );
    s0->PutScalar( 1.0 );
    Teuchos::RCP<Epetra_Vector> t0 = Teuchos::rcp( new Epetra_Vector(map) );
    jac->Apply( *s0, *t0 );
    // check norm
    double normT0;
    t0->Norm1( &normT0 );
    TEST_FLOATING_EQUALITY( normT0, controlNormT0, 1.0e-12 );
    std::cout << std::setprecision(15);
    std::cout << normT0 << std::endl;

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
    // check norm
    double normT1;
    t1->Norm1( &normT1 );
    TEST_FLOATING_EQUALITY( normT1, controlNormT1, 1.0e-12 );
    std::cout << normT1 << std::endl;

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
    // check norm
    double normT2;
    t2->Norm1( &normT2 );
    TEST_FLOATING_EQUALITY( normT2, controlNormT2, 1.0e-12 );
    std::cout << normT2 << std::endl;

    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Ginla, JacRectangleSmallHashes )
{
    std::string inputFileNameBase = "rectanglesmall";

    double mu = 1.0e-2;
    double controlNormT0 = 20.5012606103421;
    double controlNormT1 = 0.50126061034211;
    double controlNormT2 = 20.5012606103421;

    testJac( inputFileNameBase,
             mu,
             controlNormT0,
             controlNormT1,
             controlNormT2,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, JacPacmanHashes )
{
    std::string inputFileNameBase = "pacman";

    double mu = 1.0e-2;
    double controlNormT0 = 606.123970266964;
    double controlNormT1 = 0.713664749303349;
    double controlNormT2 = 605.759066191324;

    testJac( inputFileNameBase,
             mu,
             controlNormT0,
             controlNormT1,
             controlNormT2,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, JacCubeSmallHashes )
{
    std::string inputFileNameBase = "cubesmall";

    double mu = 1.0e-2;
    double controlNormT0 = 20.2913968894543;
    double controlNormT1 = 0.289990630357598;
    double controlNormT2 = 20.2899906303576;

    testJac( inputFileNameBase,
             mu,
             controlNormT0,
             controlNormT1,
             controlNormT2,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, JacCubeLargeHashes )
{
    std::string inputFileNameBase = "cubelarge";

    double mu = 1.0e-2;
    double controlNormT0 = 2006.09839315077;
    double controlNormT1 = 5.74651996372481;
    double controlNormT2 = 2005.74651996367;

    testJac( inputFileNameBase,
             mu,
             controlNormT0,
             controlNormT1,
             controlNormT2,
             out,
             success );
}
// ============================================================================
} // namespace
