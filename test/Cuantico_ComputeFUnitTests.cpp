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

#include "Cuantico_StkMesh.hpp"
#include "Cuantico_StkMeshReader.hpp"
#include "Cuantico_ScalarPotential_Constant.hpp"
#include "Cuantico_MagneticVectorPotential_ExplicitValues.hpp"
#include "Cuantico_ModelEvaluator.hpp"

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
    Cuantico::StkMeshRead( *eComm, inputFileName, data );

    // Cast the data into something more accessible.
    Teuchos::RCP<Cuantico::StkMesh> & mesh = data.get( "mesh", Teuchos::RCP<Cuantico::StkMesh>() );
    Teuchos::RCP<Epetra_Vector> & z = data.get( "psi", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::RCP<Epetra_MultiVector> & mvpValues = data.get( "A", Teuchos::RCP<Epetra_MultiVector>() );
    Teuchos::RCP<Epetra_Vector> & thickness = data.get( "thickness", Teuchos::RCP<Epetra_Vector>() );
    Teuchos::ParameterList & problemParameters = data.get( "Problem parameters", Teuchos::ParameterList() );

    // create parameter vector
    problemParameters.set("g", 1.0);
    problemParameters.set("mu", mu);

    Teuchos::RCP<Cuantico::MagneticVectorPotential::Virtual> mvp =
      Teuchos::rcp(new Cuantico::MagneticVectorPotential::ExplicitValues(mesh, mvpValues, mu));

    Teuchos::RCP<Cuantico::ScalarPotential::Virtual> sp =
      Teuchos::rcp(new Cuantico::ScalarPotential::Constant(-1.0));

    Teuchos::RCP<Cuantico::ModelEvaluator> modelEval =
        Teuchos::rcp(new Cuantico::ModelEvaluator(mesh, 1.0, sp, mvp, thickness, z));

    // Get a finite-difference approximation of df/dp.
    EpetraExt::ModelEvaluator::InArgs inArgs = modelEval->createInArgs();
    inArgs.set_x( z );
    EpetraExt::ModelEvaluator::OutArgs outArgs = modelEval->createOutArgs();

    // get parameter vector and names
    Teuchos::RCP<Epetra_Vector> p =
        Teuchos::rcp(new Epetra_Vector(*modelEval->get_p_map(0)));
    Teuchos::RCP<const Teuchos::Array<std::string> > pNames =
        modelEval->get_p_names(0);

    // create parameter vector
    for (int k=0; k<p->MyLength(); k++)
       (*p)[k] = problemParameters.get<double>( (*pNames)[k] );
    inArgs.set_p(0, p);
    Teuchos::RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(z->Map()));
    outArgs.set_f( f );
    modelEval->evalModel(inArgs, outArgs);

    // check the norms
    double normOne;
    f->Norm1( &normOne );
    TEST_FLOATING_EQUALITY( normOne, controlNormOne, 1.0e-10 );

    double normTwo;
    f->Norm2( &normTwo );
    TEST_FLOATING_EQUALITY( normTwo, controlNormTwo, 1.0e-10 );

    double normInf;
    f->NormInf( &normInf );
    TEST_FLOATING_EQUALITY( normInf, controlNormInf, 1.0e-10 );


    return;
}
// ===========================================================================
TEUCHOS_UNIT_TEST( Cuantico, ComputeFRectangleSmallHashes )
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
TEUCHOS_UNIT_TEST( Cuantico, ComputeFPacmanHashes )
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
TEUCHOS_UNIT_TEST( Cuantico, ComputeFCubeSmallHashes )
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
TEUCHOS_UNIT_TEST( Cuantico, ComputeFBrickWHoleHashes )
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
