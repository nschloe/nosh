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
  const std::map<std::string, boost::any> linear_solver_params = {
      {"package", "Belos"},
      {"method", "Pseudo Block GMRES"},
      {"parameters", list{
        {"Output Frequency", 1},
        {"Output Style", 1},
        {"Verbosity", 33}
      }}
      };

  auto init_x = std::make_shared<nosh::function>(mesh);
  init_x->putScalar(0.0);
  const auto model = std::make_shared<nosh::model>(
      mesh, init_x, f, jac, dfdp, linear_solver_params
      );


  // Check out
  // https://trilinos.org/docs/dev/packages/nox/doc/html/parameters.html
  // for a full parameter description.
  nosh::nonlinear_solve(
      model,
      {
        {"method", "Newton"},
        {"NOX", list{
          {"Status Tests", list{
            {"Test Type", "Combo"},
            {"Combo Type", "OR"},
            {"Number of Tests", 2},
            {"Test 0", list{
              {"Test Type", "Combo"},
              {"Combo Type", "AND"},
              {"Number of Tests", 2},
              {"Test 0", list{
                {"Test Type", "NormF"},
                {"Norm Type", "Two Norm"},
                {"Scale Type", "Scaled"},
                {"Tolerance", 1.0e-8}
              }},
              {"Test 1", list{
                {"Test Type", "NormWRMS"},
                {"Absolute Tolerance", 1.0e-6},
                {"Relative Tolerance", 1.0e-6}
              }},
            }},
            {"Test 1", list {
              {"Test Type", "MaxIters"},
              {"Maximum Iterations", 10}
            }}
          }},
          {"Printing", list{
           {"Output Information", list{
             {"Details", true},
             {"Outer Iteration", true},
             {"Outer Iteration Status Test", true},
             {"Inner Iteration", true},
             {"Linear Solver Details", true},
             {"Parameters", true},
             {"Warning", true},
             {"Debug", true},
             {"Test Details", true},
             {"Error", true},
             {"Stepper Iteration", true},
             {"Stepper Details", true},
             {"Stepper Parameters", true}
           }}
          }}
        }},
      }
      );

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

#endif
  nosh::write(*init_x, "out.h5m");

  return EXIT_SUCCESS;
}
