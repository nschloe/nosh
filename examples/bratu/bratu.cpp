#include "bratu.hpp"
#include <nosh.hpp>
#include <mikado.hpp>

using dict = std::map<std::string, boost::any>;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  const auto f = std::make_shared<bratu::f>(mesh);
  const auto jac = std::make_shared<bratu::jacobian>(mesh);
  const auto dfdp = std::make_shared<bratu::dfdp>(mesh);

  // const auto problem = bratu::bratu(mesh);

  // Create a model evaluator.
  const std::map<std::string, boost::any> linear_solver_params = {
      {"package", std::string("Belos")},
      {"method", std::string("Pseudo Block GMRES")},
      {"parameters", dict{
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
  const auto solution = mikado::nonlinear_solve(
      model,
      {
        {"method", std::string("Newton")},
        {"NOX", dict{
          {"Status Tests", dict{
            {"Test Type", std::string("Combo")},
            {"Combo Type", std::string("OR")},
            {"Number of Tests", 2},
            {"Test 0", dict{
              {"Test Type", std::string("Combo")},
              {"Combo Type", std::string("AND")},
              {"Number of Tests", 2},
              {"Test 0", dict{
                {"Test Type", std::string("NormF")},
                {"Norm Type", std::string("Two Norm")},
                {"Scale Type", std::string("Scaled")},
                {"Tolerance", 1.0e-8}
              }},
              {"Test 1", dict{
                {"Test Type", std::string("NormWRMS")},
                {"Absolute Tolerance", 1.0e-6},
                {"Relative Tolerance", 1.0e-6}
              }},
            }},
            {"Test 1", list {
              {"Test Type", std::string("MaxIters")},
              {"Maximum Iterations", 10}
            }}
          }},
          {"Printing", dict{
           {"Output Information", dict{
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
  // mikado::nonlinear_solve(
  //     model, x,
  //     {
  //       {"method", "Newton"},
  //       {"Direction", dict{
  //         {"method": "steepest descend"},
  //         {"forcing term alpha", 1.0}
  //       }},
  //       {"parameters", dict{
  //         {"method", "polynomial"},
  //         {"max iters", 10},
  //         {"interpolation type", "quadratic"}
  //       }}
  //     }
  //     );

  /*
  mikado::parameter_continuation(
      model, x,
      {
        {"Stepper", dict{
          {"Continuation Method", "Arc Length"},
          {"Contunuation Parameter", "mu"},
          {"Initial Value", 0.0},
          {"Max Value", 10.0},
          {"Max Nonlinear Iterations", 5},
        }};
        {"Step Size", dict{
          {"Initial Step Size", 1.0e-2}
        }},
        {"jacobian solve", dict{
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
  nosh::write(solution, mesh, "out.h5m");

  return EXIT_SUCCESS;
}
