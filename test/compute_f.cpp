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

#include <nosh.hpp>

#include <Teuchos_UnitTestHarness.hpp>

namespace
{
// =============================================================================
void
testComputeF(
    const std::string & input_filename_base,
    const double mu,
    const double control_norm_1,
    const double control_norm_2,
    const double control_norm_inf,
    Teuchos::FancyOStream & out,
    bool & success
    )
{
  // Read the data from the file.
  auto comm =  Teuchos::DefaultComm<int>::getComm();
  const std::string input_filename = "data/" + input_filename_base + ".h5m";

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
          "mu"
          ));

  // Create in_args.x
  auto in_args = model_eval->createInArgs();
  in_args.set_x(Thyra::createVector(Teuchos::rcp(z), model_eval->get_f_space()));
  // in_args.p
  auto p = Thyra::createMember(model_eval->get_p_space(0));
  Thyra::set_ele(0, 1.0, p()); // g
  Thyra::set_ele(1, 0.01, p()); // mu
  in_args.set_p(0, p);

  // Create out_args.
  auto out_args = model_eval->createOutArgs();
  auto f = Thyra::createMember(model_eval->get_f_space());
  out_args.set_f(f);

  // Fetch.
  model_eval->evalModel(in_args, out_args);

  // check the norms
  TEST_FLOATING_EQUALITY(
      Thyra::norm_1(*f),
      control_norm_1,
      1.0e-10
      );
  TEST_FLOATING_EQUALITY(
      Thyra::norm_2(*f),
      control_norm_2,
      1.0e-10
      );
  TEST_FLOATING_EQUALITY(
      Thyra::norm_inf(*f),
      control_norm_inf,
      1.0e-10
      );

  return;
}
// ===========================================================================
//TEUCHOS_UNIT_TEST(nosh, ComputeFRectangleSmallHashes)
//{
//  std::string input_filename_base = "rectanglesmall";
//
//  double mu = 1.0e-2;
//  double control_norm_1 = 0.50126061034211067;
//  double control_norm_2 = 0.24749434381636057;
//  double control_norm_inf = 0.12373710977782607;
//
//  testComputeF(input_filename_base,
//               mu,
//               control_norm_1,
//               control_norm_2,
//               control_norm_inf,
//               out,
//               success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, ComputeFPacmanHashes)
{
  std::string input_filename_base = "pacman";

  double mu = 1.0e-2;
  double control_norm_1 = 0.71366475047893463;
  double control_norm_2 = 0.12552206259336218;
  double control_norm_inf = 0.055859319123267033;

  testComputeF(input_filename_base,
               mu,
               control_norm_1,
               control_norm_2,
               control_norm_inf,
               out,
               success);
}
// ============================================================================
//TEUCHOS_UNIT_TEST(nosh, ComputeFCubeSmallHashes)
//{
//  std::string input_filename_base = "cubesmall";
//
//  double mu = 1.0e-2;
//  double control_norm_1 = 8.3541623156163313e-05;
//  double control_norm_2 = 2.9536515963905867e-05;
//  double control_norm_inf = 1.0468744547749431e-05;
//
//  testComputeF(input_filename_base,
//               mu,
//               control_norm_1,
//               control_norm_2,
//               control_norm_inf,
//               out,
//               success);
//}
// ============================================================================
TEUCHOS_UNIT_TEST(nosh, ComputeFBrickWHoleHashes)
{
  std::string input_filename_base = "brick-w-hole";

  double mu = 1.0e-2;
  double control_norm_1 = 1.8084716102419285;
  double control_norm_2 = 0.15654267585120338;
  double control_norm_inf = 0.03074423493622647;

  testComputeF(input_filename_base,
               mu,
               control_norm_1,
               control_norm_2,
               control_norm_inf,
               out,
               success);
}
// ============================================================================
} // namespace
