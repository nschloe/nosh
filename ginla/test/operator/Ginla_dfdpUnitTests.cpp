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
//#include "Ginla_KeoFactory.hpp"
#include "Ginla_ModelEvaluator.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace {

// =============================================================================
void
testDfdp( const std::string & inputFileNameBase,
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
    Teuchos::RCP<Ginla::StkMesh>     & mesh = data.get( "mesh", Teuchos::RCP<Ginla::StkMesh>() );
    Teuchos::RCP<Epetra_Vector>      & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::RCP<const Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<const Epetra_MultiVector>() );
    Teuchos::RCP<Epetra_Vector>      & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::ParameterList           & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

    // create parameter vector
    problemParameters.set( "mu", mu );
    problemParameters.set( "theta", 0.0 );
    problemParameters.set( "T", 0.0 );

    Teuchos::RCP<Ginla::MagneticVectorPotential> mvp;
    mvp = Teuchos::rcp ( new Ginla::MagneticVectorPotential ( mesh, mvpValues, problemParameters.get<double>("mu") ) );

    Teuchos::RCP<Ginla::ModelEvaluator> modelEval =
        Teuchos::rcp( new Ginla::ModelEvaluator( mesh, problemParameters, thickness, mvp, z ) );

    // Get a finite-difference approximation of df/dp.
    EpetraExt::ModelEvaluator::InArgs inArgs = modelEval->createInArgs();
    inArgs.set_x( z );
    EpetraExt::ModelEvaluator::OutArgs outArgs = modelEval->createOutArgs();

    // get paramter vector and names
    Teuchos::RCP<Epetra_Vector> p =
        Teuchos::rcp( new Epetra_Vector( *modelEval->get_p_map(0) ) );
    Teuchos::RCP<const Teuchos::Array<std::string> > pNames =
        modelEval->get_p_names(0);

    // create parameter vector
    double eps = 1.0e-6;
    problemParameters.set( "mu", mu-eps );
    for ( int k=0; k<p->MyLength(); k++ )
       (*p)[k] = problemParameters.get<double>( (*pNames)[k] );
    inArgs.set_p( 0, p );
    Teuchos::RCP<Epetra_Vector> f0 = Teuchos::rcp( new Epetra_Vector( z->Map() ) );
    outArgs.set_f( f0 );
    modelEval->evalModel( inArgs, outArgs );

    problemParameters.set( "mu", mu+eps );
    for ( int k=0; k<p->MyLength(); k++ )
       (*p)[k] = problemParameters.get<double>( (*pNames)[k] );
    inArgs.set_p( 0, p );
    Teuchos::RCP<Epetra_Vector> f1 = Teuchos::rcp( new Epetra_Vector( z->Map() ) );
    outArgs.set_f( f1 );
    modelEval->evalModel( inArgs, outArgs );

    // Calculate the finite difference approx for df/dp.
    f1->Update( -1.0, *f0, 1.0 );
    f1->Scale( 0.5/eps );

    // Get the actual derivative.
    problemParameters.set( "mu", mu );
    int muIndex = 0;
    for ( int k=0; k<p->MyLength(); k++ )
       (*p)[k] = problemParameters.get<double>( (*pNames)[k] );
    inArgs.set_p( 0, p );
    Teuchos::RCP<Epetra_Vector> dfdp = Teuchos::rcp( new Epetra_Vector( z->Map() ) );
    Teuchos::Array<int> paramIndices( Teuchos::tuple(muIndex) );
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

    double r[1];
    f1->NormInf( r );
    TEST_COMPARE( r[0], <, 1.0e-8 );

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
TEUCHOS_UNIT_TEST( Ginla, DfdpCubeLargeHtashes )
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
