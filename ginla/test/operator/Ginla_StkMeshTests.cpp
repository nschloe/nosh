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
testMesh( const std::string & inputFileNameBase,
          const unsigned int controlNumNodes,
          const double controlVolNormOne,
          const double controlVolNormTwo,
          const double controlVolNormInf,
          const Teuchos::ArrayRCP<double> & controlEdgeCoeffs,
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

    const unsigned int numNodes = mesh->getNumNodes();
//    mesh->computeFvmEntities_();
    TEST_EQUALITY( numNodes, controlNumNodes );

    const Teuchos::RCP<const Epetra_Vector> controlVols = mesh->getControlVolumes();
    double r[1];
    controlVols->Norm1( r );
    TEST_FLOATING_EQUALITY( r[0], controlVolNormOne, 1.0e-12 );
    controlVols->Norm2( r );
    TEST_FLOATING_EQUALITY( r[0], controlVolNormTwo, 1.0e-12 );
    controlVols->NormInf( r );
    TEST_FLOATING_EQUALITY( r[0], controlVolNormInf, 1.0e-12 );

    if ( eComm->NumProc() == 1 )
    {
        const Teuchos::ArrayRCP<DoubleVector> edgeCoeffs =
            mesh->getEdgeCoefficientsFallback();

        // Cannot use TEST_COMPARE_FLOATING_ARRAYS() as some
        // array values are 0.0 which messes up the relative
        // error.
        for (unsigned int k=0; k<edgeCoeffs[0].length(); k++ )
            TEST_COMPARE( fabs(edgeCoeffs[0][k]-controlEdgeCoeffs[k]), <, 1.0e-12 );
    }

    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Ginla, MeshRectangleSmallHashes )
{
    std::string inputFileNameBase = "rectanglesmall";

    unsigned int numNodes = 4;
    double controlVolNormOne = 10.0;
    double controlVolNormTwo = 5.0;
    double controlVolNormInf = 2.5;
    Teuchos::ArrayRCP<double> edgeCoeffs0(3);
    edgeCoeffs0[0] = 0.0;
    edgeCoeffs0[1] = 5.0;
    edgeCoeffs0[2] = 0.05;

    testMesh( inputFileNameBase,
              numNodes,
              controlVolNormOne,
              controlVolNormTwo,
              controlVolNormInf,
              edgeCoeffs0,
              out,
              success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, MeshPacmanHashes )
{
    std::string inputFileNameBase = "pacman";

    unsigned int numNodes = 409;
    double controlVolNormOne = 302.5227007210103;
    double controlVolNormTwo = 15.38575790933914;
    double controlVolNormInf = 1.127797467043659;
    Teuchos::ArrayRCP<double> edgeCoeffs0(3);
    edgeCoeffs0[0] = 0.1029663033064168;
    edgeCoeffs0[1] = 0.4075254593854971;
    edgeCoeffs0[2] = 0.407525890029648;

    testMesh( inputFileNameBase,
              numNodes,
              controlVolNormOne,
              controlVolNormTwo,
              controlVolNormInf,
              edgeCoeffs0,
              out,
              success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, MeshCubeSmallHashes )
{
    std::string inputFileNameBase = "cubesmall";

    unsigned int numNodes = 8;
    double controlVolNormOne = 10.0;
    double controlVolNormTwo = 3.535533905932738;
    double controlVolNormInf = 1.25;
    Teuchos::ArrayRCP<double> edgeCoeffs0(6);
    edgeCoeffs0[0] = 1.666666666666667;
    edgeCoeffs0[1] = 1.666666666666667;
    edgeCoeffs0[2] = 0.01666666666666689;
    edgeCoeffs0[3] = 0.0;
    edgeCoeffs0[4] = 0.0;
    edgeCoeffs0[5] = 0.0;

    testMesh( inputFileNameBase,
              numNodes,
              controlVolNormOne,
              controlVolNormTwo,
              controlVolNormInf,
              edgeCoeffs0,
              out,
              success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Ginla, MeshCubeLargeHashes )
{
    std::string inputFileNameBase = "cubelarge";

    unsigned int numNodes = 27000;
    double controlVolNormOne = 999.9999999999114;
    double controlVolNormTwo = 6.238401064581688;
    double controlVolNormInf = 0.04100211456823565;
    Teuchos::ArrayRCP<double> edgeCoeffs0(6);
    edgeCoeffs0[0] = 0.05747127532958986;
    edgeCoeffs0[1] = 0.05747127532958986;
    edgeCoeffs0[2] = 0.05747127532958986;
    edgeCoeffs0[3] = 0.0;
    edgeCoeffs0[4] = 0.0;
    edgeCoeffs0[5] = 0.0;

    testMesh( inputFileNameBase,
              numNodes,
              controlVolNormOne,
              controlVolNormTwo,
              controlVolNormInf,
              edgeCoeffs0,
              out,
              success );
}
// ============================================================================
} // namespace
