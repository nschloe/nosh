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

#include "Cuantico_StkMesh.hpp"
#include "Cuantico_StkMeshReader.hpp"
#include "Cuantico_MagneticVectorPotential_ExplicitValues.hpp"
#include "Cuantico_KeoContainer.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace {

// =============================================================================
void
testKeo( const std::string & inputFileNameBase,
         const double psiControlNormOne,
         const double psiControlNormInf,
         const Teuchos::Array<double> & mvpControlNormsInf,
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
    //if ( eComm->NumProc() == 1 )
      inputFileName = inputFileNameBase + ".e";
    //else
      //inputFileName = inputFileNameBase + "-balanced.par";
    // =========================================================================
    // Read the data from the file.
    Teuchos::ParameterList data;
    Cuantico::StkMeshRead( *eComm, inputFileName, data );

    // Cast the data into something more accessible.
    Teuchos::RCP<Cuantico::StkMesh> & mesh = data.get( "mesh", Teuchos::RCP<Cuantico::StkMesh>() );
    Teuchos::RCP<Epetra_Vector> & psi = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );

    // Check psi.
    double r;
    psi->Norm1( &r );
    TEST_FLOATING_EQUALITY(r, psiControlNormOne, 1.0e-12);
    psi->NormInf( &r );
    TEST_FLOATING_EQUALITY(r, psiControlNormInf, 1.0e-12);

    // Check MVP.
    // Only check the infinity-norm here as all other norms
    // only apply to vectors with non-overlapping maps.
    Teuchos::Array<double> R(3);
    TEUCHOS_ASSERT_EQUALITY(0, mvpValues->NormInf(R.getRawPtr()));
    TEST_COMPARE_FLOATING_ARRAYS(R, mvpControlNormsInf, 1.0e-12);

    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Cuantico, KeoRectangleSmallHashes )
{
    std::string inputFileNameBase = "rectanglesmall";

    const double psiControlNormOne = 4.0;
    const double psiControlNormInf = 1.0;
    Teuchos::Array<double> mvpControlNormsInf(3);
    mvpControlNormsInf[0] = 0.25;
    mvpControlNormsInf[1] = 2.5;
    mvpControlNormsInf[2] = 0.0;

    testKeo(inputFileNameBase,
            psiControlNormOne,
            psiControlNormInf,
            mvpControlNormsInf,
            out,
            success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Cuantico, KeoPacmanHashes )
{
    std::string inputFileNameBase = "pacman";

    const double psiControlNormOne = 409.0;
    const double psiControlNormInf = 1.0;
    Teuchos::Array<double> mvpControlNormsInf(3);
    mvpControlNormsInf[0] = 4.999111652374270;
    mvpControlNormsInf[1] = 5.0;
    mvpControlNormsInf[2] = 0.0;

    testKeo(inputFileNameBase,
            psiControlNormOne,
            psiControlNormInf,
            mvpControlNormsInf,
            out,
            success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Cuantico, KeoCubeSmallHashes )
{
    std::string inputFileNameBase = "cubesmall";

    const double psiControlNormOne = 8.0;
    const double psiControlNormInf = 1.0;
    Teuchos::Array<double> mvpControlNormsInf(3);
    mvpControlNormsInf[0] = 0.25;
    mvpControlNormsInf[1] = 0.25;
    mvpControlNormsInf[2] = 0.0;

    testKeo(inputFileNameBase,
            psiControlNormOne,
            psiControlNormInf,
            mvpControlNormsInf,
            out,
            success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Cuantico, KeoBrickWHoleHashes )
{
    std::string inputFileNameBase = "brick-w-hole";

    const double psiControlNormOne = 744.0;
    const double psiControlNormInf = 1.0;
    Teuchos::Array<double> mvpControlNormsInf(3);
    mvpControlNormsInf[0] = 2.5;
    mvpControlNormsInf[1] = 2.5;
    mvpControlNormsInf[2] = 0.0;

    testKeo(inputFileNameBase,
            psiControlNormOne,
            psiControlNormInf,
            mvpControlNormsInf,
            out,
            success );
}
// ============================================================================
} // namespace