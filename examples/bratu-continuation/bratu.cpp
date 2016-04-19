#include "bratu.hpp"
#include <nosh.hpp>

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
      {"package", "Belos"},
      {"method", "Pseudo Block GMRES"},
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
  // https://trilinos.org/docs/dev/packages/nox/doc/html/loca_parameters.html
  // for a full parameter description.
  nosh::parameter_continuation(
      model,
      {
        {"NOX", dict{
          {"Status Tests", dict{
            {"Test Type", "NormF"},
            {"Norm Type", "Two Norm"},
            {"Tolerance", 1.0e-8}
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
        {"LOCA", dict{
          {"Predictor", dict{
            {"Method", "Tangent"}
          }},
          {"Stepper", dict{
            {"Continuation Method", "Arc Length"},
            {"Continuation Parameter", "lmbda"},
            {"Initial Value", 2.0e-3},
            {"Min Value", -1.0},
            {"Max Value", 1.0},
            {"Max Nonlinear Iterations", 5},
          }},
          {"Step Size", dict{
            {"Initial Step Size", 1.0e-3},
            {"Min Step Size", 1.0e-5},
            {"Max Step Size", 1.0e-1},
            {"Aggressiveness", 0.1}
          }}
        }}
      }
      );

  return EXIT_SUCCESS;
}
