#include <catch.hpp>

#include <string>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>

#include <nosh.hpp>

// =============================================================================
void
testComputeF(
    const std::string & input_filename_base,
    const double mu,
    const double control_norm_1,
    const double control_norm_2,
    const double control_norm_inf
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
  auto z = mesh->get_complex_vector("psi");

  // Set the thickness field.
  auto thickness = std::make_shared<nosh::scalar_field::constant>(*mesh, 1.0);
  auto mvp = std::make_shared<nosh::vector_field::explicit_values>(*mesh, "A", mu);
  auto sp = std::make_shared<nosh::scalar_field::constant>(*mesh, -1.0);

  auto model_eval = Teuchos::rcp(new nosh::model_evaluator::nls(
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
  REQUIRE(Thyra::norm_1(*f) == Approx(control_norm_1));
  REQUIRE(Thyra::norm_2(*f) == Approx(control_norm_2));
  REQUIRE(Thyra::norm_inf(*f) == Approx(control_norm_inf));

  return;
}
// ===========================================================================
#if 0
TEST_CASE("F(x) for rectangle mesh", "[rectangle]")
{
  testComputeF(
      "rectanglesmall",
      1.0e-2,
      0.50126061034211067,
      0.24749434381636057,
      0.12373710977782607
      );
}
#endif
// ============================================================================
TEST_CASE("F(x) for pacman mesh", "[pacman]")
{
  testComputeF(
      "pacman",
      1.0e-2,
      0.71366475047893463,
      0.12552206259336218,
      0.055859319123267033
      );
}
// ============================================================================
#if 0
TEST_CASE("F(x) for cube mesh", "[cube]")
{
  testComputeF(
      "cubesmall",
      1.0e-2,
      8.3541623156163313e-05,
      2.9536515963905867e-05,
      1.0468744547749431e-05
      );
}
#endif
// ============================================================================
TEST_CASE("F(x) for brick mesh", "[brick]")
{
  testComputeF(
      "brick-w-hole",
      1.0e-2,
      1.8084716102419285,
      0.15654267585120338,
      0.03074423493622647
      );
}
// ============================================================================
