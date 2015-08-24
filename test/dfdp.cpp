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

#include <nosh.hpp>

#include <Teuchos_UnitTestHarness.hpp>

namespace
{
// ===========================================================================
void
computeFiniteDifference_(
    const Thyra::ModelEvaluator<double> & model_eval,
    const Teuchos::RCP<Thyra::VectorBase<double>> & x,
    Teuchos::RCP<Thyra::VectorBase<double> > & p,
    const int param_index,
    const Teuchos::RCP<Thyra::VectorBase<double> > & fdiff
    )
{
  const double eps = 1.0e-8;
  auto pp = p->clone_v();

  const double orig_value = Thyra::get_ele(*pp, param_index);

  auto in_args = model_eval.createInArgs();
  in_args.set_x(x);

  auto out_args = model_eval.createOutArgs();

  // Get vector at x-eps.
  Thyra::set_ele(param_index, orig_value - eps, pp());
  in_args.set_p(0, pp);
  auto f0 = Thyra::createMember(fdiff->space());
  out_args.set_f(f0);
  model_eval.evalModel(in_args, out_args);

  // Get vector at x+eps.
  Thyra::set_ele(param_index, orig_value + eps, pp());
  in_args.set_p(0, pp);
  out_args.set_f(fdiff);
  model_eval.evalModel(in_args, out_args);

  // Calculate the finite difference approx for df/dp.
  Thyra::Vp_StV(fdiff(), -1.0, *f0);
  Thyra::scale(0.5/eps, fdiff());

  return;
}
// =============================================================================
void
test_dfdp(
    const std::string & input_filename_base,
    const double mu,
    Teuchos::FancyOStream & out,
    bool & success
    )
{
  std::string input_filename = "data/" + input_filename_base + ".e";

  // Read the data from the file.
  auto mesh = nosh::read(input_filename);

  // Cast the data into something more accessible.
  auto z = mesh->get_complex_vector("psi");

  // Set the thickness field.
  auto thickness = std::make_shared<nosh::scalar_field::constant>(*mesh, 1.0);
  auto mvp = std::make_shared<nosh::vector_field::explicit_values>(*mesh, "A", mu);
  auto sp = std::make_shared<nosh::scalar_field::constant>(*mesh, -1.0);

  Teuchos::RCP<Thyra::ModelEvaluator<double>> model_eval =
    Teuchos::rcp(new nosh::model_evaluator::nls(
          mesh,
          mvp,
          sp,
          1.0,
          thickness,
          z,
          "g"
          ));

  auto vector_space_x = model_eval->get_x_space();
  auto vectorSpaceF = model_eval->get_f_space();

  // Perform the finite difference test for all parameters present in the
  // system.
  // Get a finite-difference approximation of df/dp.
  auto in_args = model_eval->createInArgs();
  auto zT = Thyra::createVector(Teuchos::rcp(z), vector_space_x);
  in_args.set_x(zT);

  auto out_args = model_eval->createOutArgs();

  // create parameter vector
  auto p = Thyra::createMember(model_eval->get_p_space(0));
  auto fdiff = Thyra::createMember(model_eval->get_f_space());

  // Get the actual derivatives.
  in_args.set_p(0, p);
  auto dfdp = Thyra::createMembers(model_eval->get_f_space(), 2);
  Thyra::ModelEvaluatorBase::Derivative<double> deriv(
      dfdp,
      Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL
      );
  out_args.set_DfDp(0, deriv);
  model_eval->evalModel(in_args, out_args);

  // Only test the first parameter "g" for now since we have to set the DKeoDP
  // above. Alternative: Use "mu" above, and take param_index 1 here.
  // TODO test both parameters
  //for (int param_index = 0; param_index < p->GlobalLength(); param_index++) {
  std::vector<int> param_indices(1);
  double r;
  for (size_t param_index = 0; param_index < 1; param_index++) {
    // Get finite difference.
    computeFiniteDifference_(
        *model_eval,
        zT,
        p,
        param_index,
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
TEUCHOS_UNIT_TEST(nosh, DfdpRectangleSmallHashes)
{
  const std::string input_filename_base = "rectanglesmall";
  const double mu = 1.0e-2;
  test_dfdp(input_filename_base, mu, out, success);
}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, DfdpPacmanHashes)
{
  const std::string input_filename_base = "pacman";
  const double mu = 1.0e-2;
  test_dfdp(input_filename_base, mu, out, success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(nosh, DfdpCubeSmallHashes)
//{
//  const std::string input_filename_base = "cubesmall";
//  const double mu = 1.0e-2;
//  test_dfdp(input_filename_base, mu, out, success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, DfdpBrickWithHoleHashes)
{
  const std::string input_filename_base = "brick-w-hole";
  const double mu = 1.0e-2;
  test_dfdp(input_filename_base, mu, out, success);
}
// ============================================================================
} // namespace
