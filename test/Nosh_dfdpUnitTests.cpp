// @HEADER
//
//    Unit test for dF/dp.
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
#include "nosh/ModelEvaluator_Nls.hpp"

#include <Teuchos_UnitTestHarness.hpp>

namespace
{
// ===========================================================================
void
computeFiniteDifference_(
    const Thyra::ModelEvaluator<double> & modelEval,
    const Teuchos::RCP<Thyra::VectorBase<double>> & x,
    Teuchos::RCP<Thyra::VectorBase<double> > & p,
    const int paramIndex,
    const Teuchos::RCP<Thyra::VectorBase<double> > & fdiff
    )
{
  const double eps = 1.0e-8;
  auto pp = p->clone_v();

  const double origValue = Thyra::get_ele(*pp, paramIndex);

  auto inArgs = modelEval.createInArgs();
  inArgs.set_x(x);

  auto outArgs = modelEval.createOutArgs();

  // Get vector at x-eps.
  Thyra::set_ele(paramIndex, origValue - eps, pp());
  inArgs.set_p(0, pp);
  auto f0 = Thyra::createMember(fdiff->space());
  outArgs.set_f(f0);
  modelEval.evalModel(inArgs, outArgs);

  // Get vector at x+eps.
  Thyra::set_ele(paramIndex, origValue + eps, pp());
  inArgs.set_p(0, pp);
  outArgs.set_f(fdiff);
  modelEval.evalModel(inArgs, outArgs);

  // Calculate the finite difference approx for df/dp.
  Thyra::Vp_StV(fdiff(), -1.0, *f0);
  Thyra::scale(0.5/eps, fdiff());

  return;
}
// =============================================================================
void
testDfdp(
    const std::string & inputFileNameBase,
    const double mu,
    Teuchos::FancyOStream & out,
    bool & success
    )
{
  auto comm = Teuchos::DefaultComm<int>::getComm();

  std::string inputFileName = "data/" + inputFileNameBase + ".e";

  // Read the data from the file.
  auto mesh = std::make_shared<Nosh::StkMesh>(
      Teuchos::get_shared_ptr(comm),
      inputFileName,
      0
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
          z,
          "g"
          ));

  auto vectorSpaceX = modelEval->get_x_space();
  auto vectorSpaceF = modelEval->get_f_space();

  // Perform the finite difference test for all parameters present in the
  // system.
  // Get a finite-difference approximation of df/dp.
  auto inArgs = modelEval->createInArgs();
  auto zT = Thyra::createVector(Teuchos::rcp(z), vectorSpaceX);
  inArgs.set_x(zT);

  auto outArgs = modelEval->createOutArgs();

  // create parameter vector
  auto p = Thyra::createMember(modelEval->get_p_space(0));
  auto fdiff = Thyra::createMember(modelEval->get_f_space());

  // Get the actual derivatives.
  inArgs.set_p(0, p);
  auto dfdp = Thyra::createMembers(modelEval->get_f_space(), 2);
  Thyra::ModelEvaluatorBase::Derivative<double> deriv(
      dfdp,
      Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL
      );
  outArgs.set_DfDp(0, deriv);
  modelEval->evalModel(inArgs, outArgs);

  // Only test the first parameter "g" for now since we have to set the DKeoDP
  // above. Alternative: Use "mu" above, and take paramIndex 1 here.
  // TODO test both parameters
  //for (int paramIndex = 0; paramIndex < p->GlobalLength(); paramIndex++) {
  std::vector<int> paramIndices(1);
  double r;
  for (size_t paramIndex = 0; paramIndex < 1; paramIndex++) {
    // Get finite difference.
    computeFiniteDifference_(
        *modelEval,
        zT,
        p,
        paramIndex,
        fdiff
        );

    // Compare the two.
    Thyra::Vp_StV(fdiff(), -1.0, *dfdp->col(0));
    r = Thyra::norm_inf(*fdiff);
    TEST_COMPARE(r, <, 1.0e-7);
  }

  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(Nosh, DfdpRectangleSmallHashes)
//{
//  const std::string inputFileNameBase = "rectanglesmall";
//  const double mu = 1.0e-2;
//  testDfdp(inputFileNameBase, mu, out, success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, DfdpPacmanHashes)
{
  const std::string inputFileNameBase = "pacman";
  const double mu = 1.0e-2;
  testDfdp(inputFileNameBase, mu, out, success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(Nosh, DfdpCubeSmallHashes)
//{
//  const std::string inputFileNameBase = "cubesmall";
//  const double mu = 1.0e-2;
//  testDfdp(inputFileNameBase, mu, out, success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(Nosh, DfdpBrickWithHoleHashes)
{
  const std::string inputFileNameBase = "brick-w-hole";
  const double mu = 1.0e-2;
  testDfdp(inputFileNameBase, mu, out, success);
}
// ============================================================================
} // namespace
