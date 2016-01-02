#include <catch.hpp>

#include <map>
#include <string>

#include <Teuchos_DefaultComm.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>

#include <nosh.hpp>

// =============================================================================
void
  testJac(
      const std::string & input_filename_base,
      const double mu,
      const double control_sum_t0,
      const double control_sum_t1,
      const double control_sum_t2
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
  REQUIRE(Thyra::dot(*s, *Js) == Approx(control_sum_t0));

  // (b) [ 1, 0, 1, 0, ... ]
  for (int k = 0; k < s->space()->dim(); k++) {
    if (k % 2 == 0) {
      Thyra::set_ele(k, 1.0, s());
    } else {
      Thyra::set_ele(k, 0.0, s());
    }
  }
  jac->apply(Thyra::NOTRANS, *s, Js(), 1.0, 0.0);
  REQUIRE(Thyra::dot(*s, *Js) == Approx(control_sum_t1));

  // (b) [ 0, 1, 0, 1, ... ]
  for (int k = 0; k < s->space()->dim(); k++) {
    if (k % 2 == 0) {
      Thyra::set_ele(k, 0.0, s());
    } else {
      Thyra::set_ele(k, 1.0, s());
    }
  }
  jac->apply(Thyra::NOTRANS, *s, Js(), 1.0, 0.0);
  REQUIRE(Thyra::dot(*s, *Js) == Approx(control_sum_t2));

  return;
}
// ============================================================================
#if 0
TEST_CASE("Jacobian for rectangle mesh", "[rectangle]")
{
  testJac(
      "rectanglesmall",
      1.0e-2,
      20.0126243424616,
      20.0063121712308,
      0.00631217123080606
      );
}
#endif
// ============================================================================
TEST_CASE("Jacobian for pacman mesh", "[pacman]")
{
  testJac(
      "pacman",
      1.0e-2,
      605.78628672795264,
      605.41584408498682,
      0.37044264296586299
      );
}
// ============================================================================
#if 0
TEST_CASE("Jacobian for cube mesh", "[cube]")
{
  testJac(
      "cubesmall",
      1.0e-2,
      20.000167083246311,
      20.000083541623155,
      8.3541623155658495e-05
      );
}
#endif
// ============================================================================
TEST_CASE("Jacobian for brick mesh", "[brick]")
{
  testJac(
      "brick-w-hole",
      1.0e-2,
      777.70784890954064,
      777.54021614941144,
      0.16763276012921419
      );
}
// ============================================================================
