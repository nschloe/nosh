#include <nosh.hpp>
#include <mikado.hpp>

#include "different-materials.hpp"

using dict = std::map<std::string, boost::any>;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  different_materials::problem p(mesh);

  nosh::function x(mesh);
  x.putScalar(0.0);

  mikado::linear_solve(
      *(p.matrix), *(p.rhs), x,
      {
        {"package", "Belos"},
        {"method", "Pseudo Block GMRES"},
        {"parameters", dict{
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
