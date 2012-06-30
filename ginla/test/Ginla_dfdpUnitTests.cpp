#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <LOCA_Parameter_Vector.H>

#include "Ginla_StkMesh.hpp"
#include "Ginla_StkMeshReader.hpp"
#include "Ginla_ScalarPotential_Constant.hpp"
#include "Ginla_MagneticVectorPotential_ExplicitValues.hpp"
#include "Ginla_ModelEvaluator.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace {

// ===========================================================================
void
computeFiniteDifference_( const Teuchos::RCP<Ginla::ModelEvaluator> & modelEval,
                          EpetraExt::ModelEvaluator::InArgs inArgs,
                          EpetraExt::ModelEvaluator::OutArgs outArgs,
                          const Teuchos::RCP<const Epetra_Vector> & p,
                          const int paramIndex,
                          const Teuchos::RCP<Epetra_Vector> & f1
                          )
{
    const double eps = 1.0e-6;
    Teuchos::RCP<Epetra_Vector> pp = Teuchos::rcp(new Epetra_Vector(*p));

    // Store the original parameter value.
    double origValue = (*p)[paramIndex];

    // Get vector at x-eps.
    (*pp)[paramIndex] = origValue - eps;
    inArgs.set_p( 0, pp );
    Teuchos::RCP<Epetra_Vector> f0 = Teuchos::rcp( new Epetra_Vector( f1->Map() ) );
    outArgs.set_f( f0 );
    modelEval->evalModel( inArgs, outArgs );

    // Get vector at x+eps.
    (*pp)[paramIndex] = origValue + eps;
    inArgs.set_p( 0, pp );
    outArgs.set_f( f1 );
    modelEval->evalModel( inArgs, outArgs );

    // Calculate the finite difference approx for df/dp.
    f1->Update( -1.0, *f0, 1.0 );
    f1->Scale( 0.5/eps );

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
    Ginla::StkMeshRead( *eComm, inputFileName, data );

    // Cast the data into something more accessible.
    Teuchos::RCP<Ginla::StkMesh> & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );
    Teuchos::RCP<Epetra_Vector> & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );
    Teuchos::RCP<Epetra_Vector> & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::ParameterList & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

    // create parameter vector
    problemParameters.set( "g", 1.0 );
    problemParameters.set( "mu", mu );

    Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp =
      Teuchos::rcp(new Ginla::MagneticVectorPotential::ExplicitValues(mesh, mvpValues, problemParameters.get<double>("mu")));

    Teuchos::RCP<Ginla::ScalarPotential::Virtual> sp =
      Teuchos::rcp(new Ginla::ScalarPotential::Constant(-1.0));

    Teuchos::RCP<Ginla::ModelEvaluator> modelEval =
        Teuchos::rcp(new Ginla::ModelEvaluator(mesh, 1.0, sp, mvp, thickness, z));

    // Get a finite-difference approximation of df/dp.
    EpetraExt::ModelEvaluator::InArgs inArgs = modelEval->createInArgs();
    inArgs.set_x( z );
    EpetraExt::ModelEvaluator::OutArgs outArgs = modelEval->createOutArgs();

    // get parameter vector and names
    Teuchos::RCP<Epetra_Vector> p =
        Teuchos::rcp( new Epetra_Vector( *modelEval->get_p_map(0) ) );
    Teuchos::RCP<const Teuchos::Array<std::string> > pNames =
        modelEval->get_p_names(0);
    for ( int k=0; k<p->MyLength(); k++ )
       (*p)[k] = problemParameters.get<double>( (*pNames)[k] );

    // -------------------------------------------------------------------------
    // Find the index of the parameter.
    std::string paramName = "mu";
    int paramIndex;
    for (int k=0; k<p->MyLength(); k++)
    {
      if ((*pNames)[k].compare(paramName) == 0)
      {
        paramIndex = k;
        break;
      }
    }
    // Get finite difference.
    Teuchos::RCP<Epetra_Vector> f1 = Teuchos::rcp(new Epetra_Vector(z->Map()));
    computeFiniteDifference_(modelEval, inArgs, outArgs, p, paramIndex, f1);

    // Get the actual derivative.
    Teuchos::RCP<Epetra_Vector> dfdp = Teuchos::rcp( new Epetra_Vector( z->Map() ) );
    inArgs.set_p( 0, p );
    Teuchos::Array<int> paramIndices( Teuchos::tuple(paramIndex) );
    EpetraExt::ModelEvaluator::DerivativeMultiVector deriv( dfdp,
                                                 EpetraExt::ModelEvaluator::DERIV_MV_BY_COL,
                                                 paramIndices
                                               );
    outArgs.set_DfDp( 0, deriv );
    Teuchos::RCP<Epetra_Vector> nullV = Teuchos::null;
    outArgs.set_f( nullV );
    modelEval->evalModel( inArgs, outArgs );

    // compare the two
    f1->Update( -1.0, *dfdp, 1.0 );
    double r;
    f1->NormInf( &r );
    TEST_COMPARE( r, <, 1.0e-8 );
    // -------------------------------------------------------------------------
    // Find the index of the parameter.
    paramName = "g";
    for (int k=0; k<p->MyLength(); k++)
    {
      if ((*pNames)[k].compare(paramName) == 0)
      {
        paramIndex = k;
        break;
      }
    }
    // Get finite difference.
    f1 = Teuchos::rcp(new Epetra_Vector(z->Map()));
    computeFiniteDifference_(modelEval, inArgs, outArgs, p, paramIndex, f1);

    // Get the actual derivative.
    dfdp = Teuchos::rcp( new Epetra_Vector( z->Map() ) );
    inArgs.set_p( 0, p );
    paramIndices = Teuchos::tuple(paramIndex);
    deriv = EpetraExt::ModelEvaluator::DerivativeMultiVector( dfdp,
                                                 EpetraExt::ModelEvaluator::DERIV_MV_BY_COL,
                                                 paramIndices
                                               );
    outArgs.set_DfDp( 0, deriv );
    outArgs.set_f( nullV );
    modelEval->evalModel( inArgs, outArgs );

    // compare the two
    f1->Update( -1.0, *dfdp, 1.0 );
    f1->NormInf( &r );
    TEST_COMPARE( r, <, 1.0e-8 );
    // -------------------------------------------------------------------------

    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Ginla, DfdpRectangleSmallHashes )
{
    std::string inputFileNameBase = "rectanglesmall";

    double mu = 1.0e-2;

    testDfdp( inputFileNameBase,
              mu,
              out,
              success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, DfdpPacmanHashes )
{
    std::string inputFileNameBase = "pacman";

    double mu = 1.0e-2;

    testDfdp( inputFileNameBase,
              mu,
              out,
              success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, DfdpCubeSmallHashes )
{
    std::string inputFileNameBase = "cubesmall";

    double mu = 1.0e-2;

    testDfdp( inputFileNameBase,
              mu,
              out,
              success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, DfdpBrickWithHoleHashes )
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
