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
#include "Ginla_MagneticVectorPotential_ExplicitValues.hpp"
#include "Ginla_KeoFactory.hpp"
#include "Ginla_JacobianOperator.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace {

// =============================================================================
void
testJac( const std::string & inputFileNameBase,
         const double mu,
         const double controlSumT0,
         const double controlSumT1,
         const double controlSumT2,
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

    Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp;
    mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential::ExplicitValues ( mesh, mvpValues, mu ) );

    Teuchos::RCP<LOCA::ParameterVector> mvpParameters =
        Teuchos::rcp( new LOCA::ParameterVector() );
    mvpParameters->addParameter( "mu", mu );

    // create a keo factory
    Teuchos::RCP<Ginla::KeoFactory> keoFactory =
        Teuchos::rcp( new Ginla::KeoFactory( mesh, thickness, mvp ) );

    // create the jacobian operator
    Teuchos::RCP<Ginla::JacobianOperator> jac =
        Teuchos::rcp( new Ginla::JacobianOperator( mesh, thickness, keoFactory, z ) );

    double sum;
    const Epetra_Map & map = jac->OperatorDomainMap();
    Epetra_Vector s(map);
    Epetra_Vector Js(map);

    // -------------------------------------------------------------------------
    // (a) [ 1, 1, 1, ... ]
    s.PutScalar( 1.0 );
    jac->Apply( s, Js );
    s.Dot( Js, &sum );
    TEST_FLOATING_EQUALITY( sum, controlSumT0, 1.0e-12 );
    // -------------------------------------------------------------------------
    // (b) [ 1, 0, 1, 0, ... ]
    double one  = 1.0;
    double zero = 0.0;
    for ( int k=0; k<map.NumMyPoints(); k++ )
    {
        if ( map.GID(k) % 2 == 0 )
            s.ReplaceMyValues( 1, &one, &k );
        else
            s.ReplaceMyValues( 1, &zero, &k );
    }
    jac->Apply( s, Js );
    s.Dot( Js, &sum );
    TEST_FLOATING_EQUALITY( sum, controlSumT1, 1.0e-12 );
    // -------------------------------------------------------------------------
    // (b) [ 0, 1, 0, 1, ... ]
    for ( int k=0; k<map.NumMyPoints(); k++ )
    {
        if ( map.GID(k) % 2 == 0 )
            s.ReplaceMyValues( 1, &zero, &k );
        else
            s.ReplaceMyValues( 1, &one, &k );
    }
    jac->Apply( s, Js );
    s.Dot( Js, &sum );
    TEST_FLOATING_EQUALITY( sum, controlSumT2, 1.0e-12 );
    // -------------------------------------------------------------------------
    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Ginla, JacRectangleSmallHashes )
{
    std::string inputFileNameBase = "rectanglesmall";

    double mu = 1.0e-2;
    double controlSumT0 = -20.0126243424616;
    double controlSumT1 = -20.0063121712308;
    double controlSumT2 = -0.00631217123080606;

    testJac( inputFileNameBase,
             mu,
             controlSumT0,
             controlSumT1,
             controlSumT2,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, JacPacmanHashes )
{
    std::string inputFileNameBase = "pacman";

    double mu = 1.0e-2;
    double controlSumT0 = -605.786286731452;
    double controlSumT1 = -605.415844086736;
    double controlSumT2 = -0.370442644715631;

    testJac( inputFileNameBase,
             mu,
             controlSumT0,
             controlSumT1,
             controlSumT2,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, JacCubeSmallHashes )
{
    std::string inputFileNameBase = "cubesmall";

    double mu = 1.0e-2;
    double controlSumT0 = -20.0084442850419;
    double controlSumT1 = -20.0042221425209;
    double controlSumT2 = -0.00422214252093753;

    testJac( inputFileNameBase,
             mu,
             controlSumT0,
             controlSumT1,
             controlSumT2,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, JacBrickWHoleHashes )
{
    std::string inputFileNameBase = "brick-w-hole";

    double mu = 1.0e-2;
    double controlSumT0 = -777.707848909512;
    double controlSumT1 = -777.540216149396;
    double controlSumT2 = -0.167632760114663;

    testJac( inputFileNameBase,
             mu,
             controlSumT0,
             controlSumT1,
             controlSumT2,
             out,
             success );
}
// ============================================================================
} // namespace
