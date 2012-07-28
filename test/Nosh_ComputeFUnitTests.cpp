// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2012  Nico Schl\"omer
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

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <LOCA_Parameter_Vector.H>

#include "Nosh_StkMesh.hpp"
#include "Nosh_StkMeshReader.hpp"
#include "Nosh_ScalarField_Constant.hpp"
#include "Nosh_MatrixBuilder_Keo.hpp"
#include "Nosh_VectorField_ExplicitValues.hpp"
#include "Nosh_ModelEvaluator.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace {

// =============================================================================
void
testComputeF( const std::string & inputFileNameBase,
              const double mu,
              const double controlNormOne,
              const double controlNormTwo,
              const double controlNormInf,
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
    Nosh::StkMeshRead( *eComm, inputFileName, 0, data );

    // Cast the data into something more accessible.
    Teuchos::RCP<Nosh::StkMesh> & mesh =
      data.get( "mesh", Teuchos::RCP<Nosh::StkMesh>() );
    Teuchos::RCP<Epetra_Vector> & z =
      data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::RCP<const Epetra_MultiVector> & mvpValues =
      data.get( "A", Teuchos::RCP<const Epetra_MultiVector>() );
    Teuchos::RCP<Epetra_Vector> & thicknessValues =
      data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );

    // Set the thickness field.
    Teuchos::RCP<Nosh::ScalarField::Virtual> thickness =
      Teuchos::rcp(new Nosh::ScalarField::Constant(1.0));

    Teuchos::RCP<Nosh::VectorField::Virtual> mvp =
      Teuchos::rcp(new Nosh::VectorField::ExplicitValues(mesh, mvpValues, mu));
    const Teuchos::RCP<Nosh::MatrixBuilder::Virtual> matrixBuilder =
      Teuchos::rcp(new Nosh::MatrixBuilder::Keo(mesh, thickness, mvp));

    Teuchos::RCP<Nosh::ScalarField::Virtual> sp =
      Teuchos::rcp(new Nosh::ScalarField::Constant(-1.0));

    Teuchos::RCP<Nosh::ModelEvaluator> modelEval =
      Teuchos::rcp(new Nosh::ModelEvaluator(mesh, matrixBuilder, sp, 1.0, thickness, z));

    // Create inArgs. Use p_init as parameters.
    EpetraExt::ModelEvaluator::InArgs inArgs = modelEval->createInArgs();
    inArgs.set_x( z );
    inArgs.set_p(0, modelEval->get_p_init(0));

    // Create outArgs.
    EpetraExt::ModelEvaluator::OutArgs outArgs = modelEval->createOutArgs();
    Teuchos::RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(z->Map()));
    outArgs.set_f( f );

    // Fetch.
    modelEval->evalModel(inArgs, outArgs);

    // check the norms
    double normOne;
    TEUCHOS_ASSERT_EQUALITY(0, f->Norm1( &normOne ));
    TEST_FLOATING_EQUALITY( normOne, controlNormOne, 1.0e-10 );

    double normTwo;
    TEUCHOS_ASSERT_EQUALITY(0, f->Norm2( &normTwo ));
    TEST_FLOATING_EQUALITY( normTwo, controlNormTwo, 1.0e-10 );

    double normInf;
    TEUCHOS_ASSERT_EQUALITY(0, f->NormInf( &normInf ));
    TEST_FLOATING_EQUALITY( normInf, controlNormInf, 1.0e-10 );

    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Nosh, ComputeFRectangleSmallHashes )
{
    std::string inputFileNameBase = "rectanglesmall";

    double mu = 1.0e-2;
    double controlNormOne = 0.50126061034211067;
    double controlNormTwo = 0.24749434381636057;
    double controlNormInf = 0.12373710977782607;

    testComputeF( inputFileNameBase,
                  mu,
                  controlNormOne,
                  controlNormTwo,
                  controlNormInf,
                  out,
                  success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Nosh, ComputeFPacmanHashes )
{
    std::string inputFileNameBase = "pacman";

    double mu = 1.0e-2;
    double controlNormOne = 0.71366475047893463;
    double controlNormTwo = 0.12552206259336218;
    double controlNormInf = 0.055859319123267033;

    testComputeF( inputFileNameBase,
                  mu,
                  controlNormOne,
                  controlNormTwo,
                  controlNormInf,
                  out,
                  success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Nosh, ComputeFCubeSmallHashes )
{
    std::string inputFileNameBase = "cubesmall";

    double mu = 1.0e-2;
    double controlNormOne = 8.3541623156163313e-05;
    double controlNormTwo = 2.9536515963905867e-05;
    double controlNormInf = 1.0468744547749431e-05;

    testComputeF( inputFileNameBase,
                  mu,
                  controlNormOne,
                  controlNormTwo,
                  controlNormInf,
                  out,
                  success );
}
// ============================================================================
TEUCHOS_UNIT_TEST( Nosh, ComputeFBrickWHoleHashes )
{
    std::string inputFileNameBase = "brick-w-hole";

    double mu = 1.0e-2;
    double controlNormOne = 1.8084716102419285;
    double controlNormTwo = 0.15654267585120338;
    double controlNormInf = 0.03074423493622647;

    testComputeF( inputFileNameBase,
                  mu,
                  controlNormOne,
                  controlNormTwo,
                  controlNormInf,
                  out,
                  success );
}
// ============================================================================
} // namespace