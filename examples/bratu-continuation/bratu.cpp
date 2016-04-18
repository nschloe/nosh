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
  // https://trilinos.org/docs/dev/packages/nox/doc/html/loca_parameters.html
  // for a full parameter description.
  nosh::parameter_continuation(
      model,
      {
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
          }}
        }},
        {"LOCA", list{
          {"Stepper", list{
            {"Continuation Method", "Arc Length"},
            {"Continuation Parameter", "lmbda"},
            {"Initial Value", 1.0e-3},
            {"Min Value", -10.0},
            {"Max Value", 10.0},
            {"Max Nonlinear Iterations", 5},
          }},
          {"Step Size", list{
            {"Initial Step Size", 1.0e-3}
          }}
        }}
      }
      );

  return EXIT_SUCCESS;
}
