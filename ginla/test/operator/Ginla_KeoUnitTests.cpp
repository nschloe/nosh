// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2011  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
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

#include <Teuchos_UnitTestHarness.hpp>

namespace {

// =============================================================================
void
testKeo( const std::string & inputFileNameBase,
         const double mu,
         const double controlNormOne,
         const double controlNormInf,
         const double controlSum,
         const double controlSumReal,
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

    Teuchos::RCP<Ginla::KeoFactory> keoFactory =
        Teuchos::rcp( new Ginla::KeoFactory( mesh, thickness, mvp ) );

    // create the kinetic energy operator
    keoFactory->updateParameters( mvpParameters );
    const Teuchos::RCP<const Epetra_CrsMatrix> keoMatrix = keoFactory->getKeo();

    // Compute matrix norms as hashes.
    // Don't check for NormFrobenius() as this one doesn't work for matrices
    // with overlapping maps.
    double normOne = keoMatrix->NormOne();
    double normInf = keoMatrix->NormInf();

    // check the values
    TEST_FLOATING_EQUALITY( normOne, controlNormOne, 1.0e-12 );
    TEST_FLOATING_EQUALITY( normInf, controlNormInf, 1.0e-12 );

    const Epetra_Map & map = keoMatrix->DomainMap();

    // Add up all the entries of the matrix.
    Epetra_Vector e( map );
    e.PutScalar( 1.0 );
    Epetra_Vector Ke( map );
    keoMatrix->Apply( e, Ke );
    double sum[1];
    e.Dot( Ke, sum );
    TEST_FLOATING_EQUALITY( sum[0], controlSum, 1.0e-10 );

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
    TEST_FLOATING_EQUALITY( sum[0], controlSumReal, 1.0e-10 );

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
    // The matrix is Hermitian, so just test that the sum of
    // the imaginary parts is (close to) 0.
    // Don't use TEST_FLOATING_EQUALITY as this one checks
    // the *relative* error.
    TEST_COMPARE( fabs(sum[0]), <, 1.0e-12 );

    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Ginla, KeoRectangleSmallHashes )
{
    std::string inputFileNameBase = "rectanglesmall";

    double mu = 1.0e-2;
    double controlNormOne = 10.2246588065616;
    double controlNormInf = 10.2246588065616;
    double controlSum     = -0.0126243424616103;
    double controlSumReal = -0.00631217123080605;

    testKeo( inputFileNameBase,
             mu,
             controlNormOne,
             controlNormInf,
             controlSum,
             controlSumReal,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, KeoPacmanHashes )
{
    std::string inputFileNameBase = "pacman";

    double mu = 1.0e-2;
    double controlNormOne = 10.0005208574565;
    double controlNormInf = 10.0005208574565;
    double controlSum     = -0.740885289431222;
    double controlSumReal = -0.370442644715617;

    testKeo( inputFileNameBase,
             mu,
             controlNormOne,
             controlNormInf,
             controlSum,
             controlSumReal,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, KeoCubeSmallHashes )
{
    std::string inputFileNameBase = "cubesmall";

    double mu = 1.0e-2;
    double controlNormOne = 10.1456474156918;
    double controlNormInf = 10.1456474156918;
    double controlSum     = -0.00844428504187249;
    double controlSumReal = -0.0042221425209367988;

    testKeo( inputFileNameBase,
             mu,
             controlNormOne,
             controlNormInf,
             controlSum,
             controlSumReal,
             out,
             success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, KeoCubeLargeHashes )
{
    std::string inputFileNameBase = "cubelarge";

    double mu = 1.0e-2;
    double controlNormOne = 4.15276626198951;
    double controlNormInf = 4.15276626198951;
    double controlSum     = -0.834815525936051;
    double controlSumReal = -0.417407762968018;

    testKeo( inputFileNameBase,
             mu,
             controlNormOne,
             controlNormInf,
             controlSum,
             controlSumReal,
             out,
             success );
}
// ============================================================================
} // namespace