#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <LOCA_Parameter_Vector.H>

#include "Cuantico_StkMesh.hpp"
#include "Cuantico_StkMeshReader.hpp"
#include "Cuantico_ScalarPotential_Constant.hpp"
#include "Cuantico_MagneticVectorPotential_ExplicitValues.hpp"
#include "Cuantico_ModelEvaluator.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace {

// ===========================================================================
void
computeFiniteDifference_( const Teuchos::RCP<Cuantico::ModelEvaluator> & modelEval,
                          EpetraExt::ModelEvaluator::InArgs inArgs,
                          EpetraExt::ModelEvaluator::OutArgs outArgs,
                          const Teuchos::RCP<const Epetra_Vector> & p,
                          const int paramIndex,
                          const Teuchos::RCP<Epetra_Vector> & fdiff
                          )
{
    const double eps = 1.0e-6;
    Teuchos::RCP<Epetra_Vector> pp = Teuchos::rcp(new Epetra_Vector(*p));

    // Store the original parameter value.
    const double origValue = (*p)[paramIndex];

    // Get vector at x-eps.
    (*pp)[paramIndex] = origValue - eps;
    inArgs.set_p( 0, pp );
    Teuchos::RCP<Epetra_Vector> f0 = Teuchos::rcp( new Epetra_Vector( fdiff->Map() ) );
    outArgs.set_f( f0 );
    modelEval->evalModel( inArgs, outArgs );

    // Get vector at x+eps.
    (*pp)[paramIndex] = origValue + eps;
    inArgs.set_p( 0, pp );
    outArgs.set_f( fdiff );
    modelEval->evalModel( inArgs, outArgs );

    // Calculate the finite difference approx for df/dp.
    TEUCHOS_ASSERT_EQUALITY(0, fdiff->Update( -1.0, *f0, 1.0 ));
    TEUCHOS_ASSERT_EQUALITY(0, fdiff->Scale( 0.5/eps ));

    return;
}
// =============================================================================
void
testDfdp(const std::string & inputFileNameBase,
         const double mu,
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
    Teuchos::RCP<Epetra_Vector> & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );
    Teuchos::RCP<Epetra_Vector> & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );

    Teuchos::RCP<Cuantico::MagneticVectorPotential::Virtual> mvp =
      Teuchos::rcp(new Cuantico::MagneticVectorPotential::ExplicitValues(mesh, mvpValues, mu));

    Teuchos::RCP<Cuantico::ScalarPotential::Virtual> sp =
      Teuchos::rcp(new Cuantico::ScalarPotential::Constant(-1.0));

    Teuchos::RCP<Cuantico::ModelEvaluator> modelEval =
      Teuchos::rcp(new Cuantico::ModelEvaluator(mesh, 1.0, sp, mvp, thickness, z));

    // -------------------------------------------------------------------------
    // Perform the finite difference test for all parameters present in the
    // system.
    // Get a finite-difference approximation of df/dp.
    EpetraExt::ModelEvaluator::InArgs inArgs = modelEval->createInArgs();
    inArgs.set_x( z );
    EpetraExt::ModelEvaluator::OutArgs outArgs = modelEval->createOutArgs();

    // Get a the initial parameter vector.
    const Teuchos::RCP<const Epetra_Vector> p = modelEval->get_p_init(0);
    Teuchos::RCP<Epetra_Vector> fdiff = Teuchos::rcp(new Epetra_Vector(z->Map()));
    Teuchos::RCP<Epetra_Vector> dfdp = Teuchos::rcp(new Epetra_Vector(z->Map()));
    EpetraExt::ModelEvaluator::DerivativeMultiVector deriv;
    Teuchos::Array<int> paramIndices(1);
    const Teuchos::RCP<Epetra_Vector> nullV = Teuchos::null;
    double r;
    for (int paramIndex=0; paramIndex < p->GlobalLength(); paramIndex++)
    {
      // Get finite difference.
      computeFiniteDifference_(modelEval, inArgs, outArgs, p, paramIndex, fdiff);

      // Get the actual derivative.
      inArgs.set_p( 0, p );
      paramIndices[0] = paramIndex;
      deriv = EpetraExt::ModelEvaluator::DerivativeMultiVector(dfdp,
                                                               EpetraExt::ModelEvaluator::DERIV_MV_BY_COL,
                                                               paramIndices);
      outArgs.set_DfDp( 0, deriv );
      outArgs.set_f( nullV );
      modelEval->evalModel( inArgs, outArgs );

      // Compare the two.
      TEUCHOS_ASSERT_EQUALITY(0, fdiff->Update( -1.0, *dfdp, 1.0 ));
      TEUCHOS_ASSERT_EQUALITY(0, fdiff->NormInf( &r ));
      TEST_COMPARE( r, <, 1.0e-8 );
    }
    // -------------------------------------------------------------------------

    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Cuantico, DfdpRectangleSmallHashes )
{
    std::string inputFileNameBase = "rectanglesmall";

    double mu = 1.0e-2;

    testDfdp( inputFileNameBase,
              mu,
              out,
              success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Cuantico, DfdpPacmanHashes )
{
    std::string inputFileNameBase = "pacman";

    double mu = 1.0e-2;

    testDfdp( inputFileNameBase,
              mu,
              out,
              success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Cuantico, DfdpCubeSmallHashes )
{
    std::string inputFileNameBase = "cubesmall";

    double mu = 1.0e-2;

    testDfdp( inputFileNameBase,
              mu,
              out,
              success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Cuantico, DfdpBrickWithHoleHashes )
{
    std::string inputFileNameBase = "brick-w-hole";

    double mu = 1.0e-2;

    testDfdp( inputFileNameBase,
              mu,
              out,
              success );
}
// ============================================================================
} // namespace
