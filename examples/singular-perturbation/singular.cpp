#include "singular.hpp"
#include <nosh.hpp>
#include <memory>

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.e");

  const auto bc1 = std::make_shared<singular::bc1>();

  singular::singular matrix(mesh, {bc1});

  const nosh::constant rhs(1.0);

  nosh::function x(mesh);

  nosh::scaled_linear_solve(
    matrix, rhs, x,
    {
      {"method", "Pseudo Block CG"},
      {
        "parameters", list{
        {"Convergence Tolerance", 1.0e-10},
        {"Output Frequency", 1},
        {"Output Style", 1},
        {"Verbosity", 33}
        }
      }
    }
  );

  nosh::write(x, "out.e");

  return EXIT_SUCCESS;
}
