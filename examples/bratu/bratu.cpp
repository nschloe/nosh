#include "bratu.hpp"
#include <nosh.hpp>

using list = std::map<std::string, boost::any>;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  const auto f = bratu::f(mesh);
  const auto jac = bratu::jacobian(mesh);

#if 0
  const auto problem = bratu::bratu(mesh);

  // Create a model evaluator.
  // Can be used at everything in Trilinos that accepts such a thing.
  const auto model = nosh::model(
      problem.f, problem.jac, problem.dfdp,
      {
        {"linear solver package", "Belos"},
        {"method", "Pseudo Block CG"},
        {"preconditioner", problem.prec}
      }
      );

  // starting value
  nosh::function x(mesh);
  x.putScalar(0.0);

  nosh::nonlinear_solve(
      model, x,
      {
        {"method", "Newton"},
      }
      );

  // nosh::nonlinear_solve(
  //     model, x,
  //     {
  //       {"method", "Newton"},
  //       {"Direction", list{
  //         {"method": "steepest descend"},
  //         {"forcing term alpha", 1.0}
  //       }},
  //       {"parameters", list{
  //         {"method", "polynomial"},
  //         {"max iters", 10},
  //         {"interpolation type", "quadratic"}
  //       }}
  //     }
  //     );

  /*
  nosh::parameter_continuation(
      model, x,
      {
        {"Stepper", list{
          {"Continuation Method", "Arc Length"},
          {"Contunuation Parameter", "mu"},
          {"Initial Value", 0.0},
          {"Max Value", 10.0},
          {"Max Nonlinear Iterations", 5},
        }};
        {"Step Size", list{
          {"Initial Step Size", 1.0e-2}
        }},
        {"jacobian solve", list{
          {"jacobian", jac},
          {"linear solver package", "Belos"},
          {"method", "Pseudo Block CG"},
          {"preconditioner", M}
        }},
        {"fpdp", dfdp}
      }
      );
      */

  nosh::write(x, "out.h5m");
#endif

  return EXIT_SUCCESS;
}
