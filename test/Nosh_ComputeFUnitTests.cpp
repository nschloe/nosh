// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2012  Nico Schlömer
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
#include <string>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>

#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Constant.hpp"
#include "nosh/ParameterMatrix_Keo.hpp"
#include "nosh/ParameterMatrix_DKeoDP.hpp"
#include "nosh/VectorField_ExplicitValues.hpp"
#include "nosh/ModelEvaluator_Nls.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace
{
// =============================================================================
void
testComputeF(
    const std::string & inputFileNameBase,
    const double mu,
    const double controlNormOne,
    const double controlNormTwo,
    const double controlNormInf,
    Teuchos::FancyOStream & out,
    bool & success
    )
{
  auto comm = Teuchos::DefaultComm<int>::getComm();

  std::string inputFileName = "data/" + inputFileNameBase + ".e";

  // Read the data from the file.
  auto mesh = std::make_shared<Nosh::StkMesh>(
      Teuchos::get_shared_ptr(comm), inputFileName, 0
      );

  // Cast the data into something more accessible.
  auto z = mesh->createComplexVector("psi");

  // Set the thickness field.
  auto thickness = std::make_shared<Nosh::ScalarField::Constant>(*mesh, 1.0);
  auto mvp = std::make_shared<Nosh::VectorField::ExplicitValues>(*mesh, "A", mu);

  auto sp = std::make_shared<Nosh::ScalarField::Constant>(*mesh, -1.0);

  Teuchos::RCP<Thyra::ModelEvaluator<double>> modelEval =
    Teuchos::rcp(new Nosh::ModelEvaluator::Nls(
          mesh,
          mvp,
          sp,
          1.0,
          thickness,
          z
          ));

  // Create inArgs.x
  auto inArgs = modelEval->createInArgs();
  inArgs.set_x(Thyra::createVector(Teuchos::rcp(z), modelEval->get_f_space()));
  // inArgs.p
  auto p = Thyra::createMember(modelEval->get_p_space(0));
  Thyra::set_ele(0, 1.0, p()); // g
  Thyra::set_ele(1, 0.01, p()); // mu
  inArgs.set_p(0, p);

  auto mymap = mesh->getComplexNonOverlapMap();
  auto space = Thyra::createVectorSpace<int,int>(Teuchos::rcp(mymap));

  // Create outArgs.
  auto outArgs = modelEval->createOutArgs();
  auto f = Thyra::createMember(modelEval->get_f_space());
  outArgs.set_f(f);

  // Fetch.
  modelEval->evalModel(inArgs, outArgs);

  // check the norms
  TEST_FLOATING_EQUALITY(
      Thyra::norm_1(*f),
      controlNormOne,
      1.0e-10
      );
  TEST_FLOATING_EQUALITY(
      Thyra::norm_2(*f),
      controlNormTwo,
      1.0e-10
      );
  TEST_FLOATING_EQUALITY(
      Thyra::norm_inf(*f),
      controlNormInf,
      1.0e-10
      );

  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(Nosh, ComputeFRectangleSmallHashes)
//{
//  std::string inputFileNameBase = "rectanglesmall";
//
//  double mu = 1.0e-2;
//  double controlNormOne = 0.50126061034211067;
//  double controlNormTwo = 0.24749434381636057;
//  double controlNormInf = 0.12373710977782607;
//
//  testComputeF(inputFileNameBase,
//               mu,
//               controlNormOne,
//               controlNormTwo,
//               controlNormInf,
//               out,
//               success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, ComputeFPacmanHashes)
{
  std::string inputFileNameBase = "pacman";

  double mu = 1.0e-2;
  double controlNormOne = 0.71366475047893463;
  double controlNormTwo = 0.12552206259336218;
  double controlNormInf = 0.055859319123267033;

  testComputeF(inputFileNameBase,
               mu,
               controlNormOne,
               controlNormTwo,
               controlNormInf,
               out,
               success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(Nosh, ComputeFCubeSmallHashes)
//{
//  std::string inputFileNameBase = "cubesmall";
//
//  double mu = 1.0e-2;
//  double controlNormOne = 8.3541623156163313e-05;
//  double controlNormTwo = 2.9536515963905867e-05;
//  double controlNormInf = 1.0468744547749431e-05;
//
//  testComputeF(inputFileNameBase,
//               mu,
//               controlNormOne,
//               controlNormTwo,
//               controlNormInf,
//               out,
//               success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, ComputeFBrickWHoleHashes)
{
  std::string inputFileNameBase = "brick-w-hole";

  double mu = 1.0e-2;
  double controlNormOne = 1.8084716102419285;
  double controlNormTwo = 0.15654267585120338;
  double controlNormInf = 0.03074423493622647;

  testComputeF(inputFileNameBase,
               mu,
               controlNormOne,
               controlNormTwo,
               controlNormInf,
               out,
               success);
}
// ============================================================================
} // namespace
