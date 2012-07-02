#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Vector.h>
#include <LOCA_Parameter_Vector.H>

#include "Cuantico_StkMesh.hpp"
#include "Cuantico_StkMeshReader.hpp"
#include "Cuantico_MagneticVectorPotential_ExplicitValues.hpp"
#include "Cuantico_ScalarPotential_Constant.hpp"
#include "Cuantico_KeoContainer.hpp"
#include "Cuantico_JacobianOperator.hpp"

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
    Cuantico::StkMeshRead( *eComm, inputFileName, data );

    // Cast the data into something more accessible.
    Teuchos::RCP<Cuantico::StkMesh> & mesh = data.get( "mesh", Teuchos::RCP<Cuantico::StkMesh>() );
    Teuchos::RCP<Epetra_Vector> & psi = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );
    Teuchos::RCP<Epetra_Vector> & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::ParameterList & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

    const double g = 1.0;
    Teuchos::Array<double> mvpParameters(1);
    mvpParameters[0] = mu;
    Teuchos::Array<double> spParameters(1);
    spParameters[0] = 0.0; // T

    Teuchos::RCP<Cuantico::MagneticVectorPotential::Virtual> mvp =
      Teuchos::rcp(new Cuantico::MagneticVectorPotential::ExplicitValues(mesh, mvpValues, mu));

    Teuchos::RCP<Cuantico::ScalarPotential::Virtual> sp =
      Teuchos::rcp(new Cuantico::ScalarPotential::Constant(-1.0));

    // create a keo factory
    Teuchos::RCP<Cuantico::KeoContainer> keoContainer =
        Teuchos::rcp(new Cuantico::KeoContainer(mesh, thickness, mvp));

    // create the jacobian operator
    Teuchos::RCP<Cuantico::JacobianOperator> jac =
      Teuchos::rcp(new Cuantico::JacobianOperator(mesh, sp, thickness, keoContainer));
    jac->rebuild(g, spParameters, mvpParameters, psi);

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
    TEST_EQUALITY(0, jac->Apply( s, Js ));
    TEST_EQUALITY(0, s.Dot( Js, &sum ));
    TEST_FLOATING_EQUALITY( sum, controlSumT2, 1.0e-10 );
    // -------------------------------------------------------------------------
    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Cuantico, JacRectangleSmallHashes )
{
    std::string inputFileNameBase = "rectanglesmall";

    double mu = 1.0e-2;
    double controlSumT0 = 20.0126243424616;
    double controlSumT1 = 20.0063121712308;
    double controlSumT2 = 0.00631217123080606;

    testJac( inputFileNameBase,
             mu,
             controlSumT0,
             controlSumT1,
             controlSumT2,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Cuantico, JacPacmanHashes )
{
    std::string inputFileNameBase = "pacman";

    double mu = 1.0e-2;
    double controlSumT0 = 605.78628672795264;
    double controlSumT1 = 605.41584408498682;
    double controlSumT2 = 0.37044264296586299;

    testJac( inputFileNameBase,
             mu,
             controlSumT0,
             controlSumT1,
             controlSumT2,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Cuantico, JacCubeSmallHashes )
{
    std::string inputFileNameBase = "cubesmall";

    double mu = 1.0e-2;
    double controlSumT0 = 20.000167083246311;
    double controlSumT1 = 20.000083541623155;
    double controlSumT2 = 8.3541623155658495e-05;

    testJac( inputFileNameBase,
             mu,
             controlSumT0,
             controlSumT1,
             controlSumT2,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Cuantico, JacBrickWHoleHashes )
{
    std::string inputFileNameBase = "brick-w-hole";

    double mu = 1.0e-2;
    double controlSumT0 = 777.70784890954064;
    double controlSumT1 = 777.54021614941144;
    double controlSumT2 = 0.16763276012921419;

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
