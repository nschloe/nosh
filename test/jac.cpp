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

#include <nosh.hpp>

#include <Teuchos_UnitTestHarness.hpp>

namespace
{

// =============================================================================
void
  testJac(
      const std::string & input_filename_base,
      const double mu,
      const double control_sum_t0,
      const double control_sum_t1,
      const double control_sum_t2,
      Teuchos::FancyOStream & out,
      bool & success
      )
{
  // Read the data from the file.
  auto comm =  Teuchos::DefaultComm<int>::getComm();
  const int size = comm->getSize();
  const std::string input_filename = (size == 1) ?
    "data/" + input_filename_base + ".h5m" :
    "data/" + input_filename_base + "-" + std::to_string(size) + ".h5m"
    ;

  // Read the data from the file.
  auto mesh = nosh::read(input_filename);

  // Cast the data into something more accessible.
  auto psi = mesh->get_complex_vector("psi");

  std::map<std::string, double> params;
  params["g"] = 1.0;
  params["mu"] = mu;

  auto mvp = std::make_shared<nosh::vector_field::explicit_values>(*mesh, "A", mu);
  auto sp = std::make_shared<nosh::scalar_field::constant>(*mesh, -1.0);
  auto thickness = std::make_shared<nosh::scalar_field::constant>(*mesh, 1.0);

  Teuchos::RCP<nosh::model_evaluator::nls> model_eval =
    Teuchos::rcp(new nosh::model_evaluator::nls(
          mesh,
          mvp,
          sp,
          1.0,
          thickness,
          psi,
          "mu"
          ));

  // set parameters
  auto in_args = model_eval->createInArgs();
  auto p = Thyra::createMember(model_eval->get_p_space(0));
  auto p_names = model_eval->get_p_names(0);
  for (int i=0; i<p_names->size(); i++) {
    Thyra::set_ele(i, params.at((*p_names)[i]), p());
  }
  in_args.set_p(0, p);
  in_args.set_x(Thyra::createVector(
        Teuchos::rcp(psi),
        model_eval->get_x_space()
        ));

  // get the jacobian from the model evaluator
  auto jac = model_eval->create_W_op();

  auto out_args = model_eval->createOutArgs();
  out_args.set_W_op(jac);

  // call the model
  model_eval->evalModel(in_args, out_args);

  TEUCHOS_ASSERT(!jac.is_null());

  auto domain = jac->domain();
  TEUCHOS_ASSERT(!domain.is_null());
  auto s = Thyra::createMember(domain);
  auto range = jac->range();
  auto Js = Thyra::createMember(range);

  // (a) [ 1, 1, 1, ... ]
  Thyra::put_scalar<double>(1.0, s());
  jac->apply(Thyra::NOTRANS, *s, Js(), 1.0, 0.0);
  TEST_FLOATING_EQUALITY(
      Thyra::dot(*s, *Js),
      control_sum_t0,
      1.0e-12
      );

  // (b) [ 1, 0, 1, 0, ... ]
  for (int k = 0; k < s->space()->dim(); k++) {
    if (k % 2 == 0) {
      Thyra::set_ele(k, 1.0, s());
    } else {
      Thyra::set_ele(k, 0.0, s());
    }
  }
  jac->apply(Thyra::NOTRANS, *s, Js(), 1.0, 0.0);
  TEST_FLOATING_EQUALITY(
      Thyra::dot(*s, *Js),
      control_sum_t1,
      1.0e-12
      );

  // (b) [ 0, 1, 0, 1, ... ]
  for (int k = 0; k < s->space()->dim(); k++) {
    if (k % 2 == 0) {
      Thyra::set_ele(k, 0.0, s());
    } else {
      Thyra::set_ele(k, 1.0, s());
    }
  }
  jac->apply(Thyra::NOTRANS, *s, Js(), 1.0, 0.0);
  TEST_FLOATING_EQUALITY(
      Thyra::dot(*s, *Js),
      control_sum_t2,
      1.0e-10
      );

  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(nosh, JacRectangleSmallHashes)
//{
//  std::string input_filename_base = "rectanglesmall";
//
//  double mu = 1.0e-2;
//  double control_sum_t0 = 20.0126243424616;
//  double control_sum_t1 = 20.0063121712308;
//  double control_sum_t2 = 0.00631217123080606;
//
//  testJac(input_filename_base,
//          mu,
//          control_sum_t0,
//          control_sum_t1,
//          control_sum_t2,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, JacPacmanHashes)
{
  std::string input_filename_base = "pacman";

  double mu = 1.0e-2;
  double control_sum_t0 = 605.78628672795264;
  double control_sum_t1 = 605.41584408498682;
  double control_sum_t2 = 0.37044264296586299;

  testJac(input_filename_base,
          mu,
          control_sum_t0,
          control_sum_t1,
          control_sum_t2,
          out,
          success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(nosh, JacCubeSmallHashes)
//{
//  std::string input_filename_base = "cubesmall";
//
//  double mu = 1.0e-2;
//  double control_sum_t0 = 20.000167083246311;
//  double control_sum_t1 = 20.000083541623155;
//  double control_sum_t2 = 8.3541623155658495e-05;
//
//  testJac(input_filename_base,
//          mu,
//          control_sum_t0,
//          control_sum_t1,
//          control_sum_t2,
//          out,
//          success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, JacBrickWHoleHashes)
{
  std::string input_filename_base = "brick-w-hole";

  double mu = 1.0e-2;
  double control_sum_t0 = 777.70784890954064;
  double control_sum_t1 = 777.54021614941144;
  double control_sum_t2 = 0.16763276012921419;

  testJac(input_filename_base,
          mu,
          control_sum_t0,
          control_sum_t1,
          control_sum_t2,
          out,
          success);
}
// ============================================================================
} // namespace
