// @HEADER
//
//    Unit tests for the Jacobian operator.
//    Copyright (C) 2012--2014  Nico Schl\"omer
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
#include <map>
#include <string>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>

#include <Thyra_LinearSolverBuilderBase.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>

#include "nosh/StkMesh.hpp"
#include "nosh/VectorField_ExplicitValues.hpp"
#include "nosh/ScalarField_Constant.hpp"
#include "nosh/ParameterMatrix_Keo.hpp"
#include "nosh/ParameterMatrix_DKeoDP.hpp"
#include "nosh/JacobianOperator.hpp"
#include "nosh/ModelEvaluator_Nls.hpp"
#include "nosh/ModelEvaluator_Nls.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace
{

// =============================================================================
void
  testJac(
      const std::string & inputFileNameBase,
      const double mu,
      const double controlSumT0,
      const double controlSumT1,
      const double controlSumT2,
      Teuchos::FancyOStream & out,
      bool & success
      )
{
  Teuchos::RCP<const Teuchos::Comm<int>> comm =
    Teuchos::DefaultComm<int>::getComm();

  std::string inputFileName = "data/" + inputFileNameBase + ".e";
  // =========================================================================
  // Read the data from the file.
  std::shared_ptr<Nosh::StkMesh> mesh(
      new Nosh::StkMesh(Teuchos::get_shared_ptr(comm), inputFileName, 0)
      );

  // Cast the data into something more accessible.
  auto psi = mesh->createComplexVector("psi");

  std::map<std::string, double> params;
  params["g"] = 1.0;
  params["mu"] = mu;

  std::shared_ptr<Nosh::VectorField::Virtual> mvp(
      new Nosh::VectorField::ExplicitValues(*mesh, "A", mu)
      );

  std::shared_ptr<Nosh::ScalarField::Virtual> sp(
      new Nosh::ScalarField::Constant(*mesh, -1.0)
      );

  // Set the thickness field.
  std::shared_ptr<Nosh::ScalarField::Virtual> thickness(
      new Nosh::ScalarField::Constant(*mesh, 1.0)
      );

  Teuchos::RCP<Nosh::ModelEvaluator::Nls> modelEval =
    Teuchos::rcp(new Nosh::ModelEvaluator::Nls(
          mesh,
          mvp,
          sp,
          1.0,
          thickness,
          psi
          ));

  // set parameters
  auto inArgs = modelEval->createInArgs();
  auto p = Thyra::createMember(modelEval->get_p_space(0));
  auto pNames = modelEval->get_p_names(0);
  for (int i=0; i<pNames->size(); i++) {
    Thyra::set_ele(i, params.at((*pNames)[i]), p());
  }
  inArgs.set_p(0, p);
  inArgs.set_x(Thyra::createVector(
        Teuchos::rcp(psi),
        modelEval->get_x_space()
        ));

  // get the jacobian from the model evaluator
  auto jac = modelEval->create_W_op();

  auto outArgs = modelEval->createOutArgs();
  outArgs.set_W_op(jac);

  // call the model
  modelEval->evalModel(inArgs, outArgs);

  TEUCHOS_ASSERT(!jac.is_null());

  auto domain = jac->domain();
  TEUCHOS_ASSERT(!domain.is_null());
  auto s = Thyra::createMember(domain);
  auto range = jac->range();
  auto Js = Thyra::createMember(range);

  double sum;
  // -------------------------------------------------------------------------
  // (a) [ 1, 1, 1, ... ]
  Thyra::put_scalar<double>(1.0, s());
  jac->apply(Thyra::NOTRANS, *s, Js(), 1.0, 0.0);
  sum = Thyra::dot(*s, *Js);
  TEST_FLOATING_EQUALITY(sum, controlSumT0, 1.0e-12);
  // -------------------------------------------------------------------------
  // (b) [ 1, 0, 1, 0, ... ]
  for (int k = 0; k < s->space()->dim(); k++) {
    if (k % 2 == 0)
      Thyra::set_ele(k, 1.0, s());
    else
      Thyra::set_ele(k, 0.0, s());
  }
  jac->apply(Thyra::NOTRANS, *s, Js(), 1.0, 0.0);
  sum = Thyra::dot(*s, *Js);
  TEST_FLOATING_EQUALITY(sum, controlSumT1, 1.0e-12);
  // -------------------------------------------------------------------------
  // (b) [ 0, 1, 0, 1, ... ]
  for (int k = 0; k < s->space()->dim(); k++) {
    if (k % 2 == 0)
      Thyra::set_ele(k, 0.0, s());
    else
      Thyra::set_ele(k, 1.0, s());
  }
  jac->apply(Thyra::NOTRANS, *s, Js(), 1.0, 0.0);
  sum = Thyra::dot(*s, *Js);
  TEST_FLOATING_EQUALITY(sum, controlSumT2, 1.0e-10);
  // -------------------------------------------------------------------------
  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(Nosh, JacRectangleSmallHashes)
//{
//  std::string inputFileNameBase = "rectanglesmall";
//
//  double mu = 1.0e-2;
//  double controlSumT0 = 20.0126243424616;
//  double controlSumT1 = 20.0063121712308;
//  double controlSumT2 = 0.00631217123080606;
//
//  testJac(inputFileNameBase,
//          mu,
//          controlSumT0,
//          controlSumT1,
//          controlSumT2,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, JacPacmanHashes)
{
  std::string inputFileNameBase = "pacman";

  double mu = 1.0e-2;
  double controlSumT0 = 605.78628672795264;
  double controlSumT1 = 605.41584408498682;
  double controlSumT2 = 0.37044264296586299;

  testJac(inputFileNameBase,
          mu,
          controlSumT0,
          controlSumT1,
          controlSumT2,
          out,
          success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(Nosh, JacCubeSmallHashes)
//{
//  std::string inputFileNameBase = "cubesmall";
//
//  double mu = 1.0e-2;
//  double controlSumT0 = 20.000167083246311;
//  double controlSumT1 = 20.000083541623155;
//  double controlSumT2 = 8.3541623155658495e-05;
//
//  testJac(inputFileNameBase,
//          mu,
//          controlSumT0,
//          controlSumT1,
//          controlSumT2,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, JacBrickWHoleHashes)
{
  std::string inputFileNameBase = "brick-w-hole";

  double mu = 1.0e-2;
  double controlSumT0 = 777.70784890954064;
  double controlSumT1 = 777.54021614941144;
  double controlSumT2 = 0.16763276012921419;

  testJac(inputFileNameBase,
          mu,
          controlSumT0,
          controlSumT1,
          controlSumT2,
          out,
          success);
}
// ============================================================================
} // namespace
