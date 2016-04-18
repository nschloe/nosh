#include <nosh.hpp>

#include "different-materials.hpp"

using list = std::map<std::string, boost::any>;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  different_materials::problem p(mesh);

  nosh::function x(mesh);
  x.putScalar(0.0);

  //nosh::scaled_linear_solve(
  nosh::linear_solve(
      p, x,
      {
        {"package", "Belos"},
        {"method", "Pseudo Block GMRES"},
        {"parameters", list{
          {"Convergence Tolerance", 1.0e-10},
          {"Output Frequency", 1},
          {"Output Style", 1},
          {"Verbosity", 33}
        }},
        {"preconditioner", "MueLu"}
      }
      );

  nosh::write(x, "out.h5m");

  return EXIT_SUCCESS;
}
