#include "bratu.hpp"
#include <nosh.hpp>

using list = std::map<std::string, boost::any>;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  const auto f = std::make_shared<bratu::f>(mesh);
  const auto jac = std::make_shared<bratu::jacobian>(mesh);
  const auto dfdp = std::make_shared<bratu::dfdp>(mesh);

  // const auto problem = bratu::bratu(mesh);

  // Create a model evaluator.
  // Can be used with everything in Trilinos that accepts such a thing.
  const auto model = std::make_shared<nosh::model>(
      mesh, f, jac, dfdp
      // {
      //   {"linear solver package", "Belos"},
      //   {"method", "Pseudo Block CG"},
      //   // {"preconditioner", problem.prec}
      // }
  );


  nosh::nonlinear_solve(
      model,
      {
        {"method", "Newton"},
      }
      );

  // // starting value
  // nosh::function x(mesh);
  // x.putScalar(0.0);

  // nosh::nonlinear_solve(
  //     model, x,
  //     {
  //       {"method", "Newton"},
  //     }
  //     );

#if 0
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
