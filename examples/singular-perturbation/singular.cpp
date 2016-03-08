#include "singular.hpp"
#include <nosh.hpp>
#include <memory>

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  const singular::singular matrix(mesh);

  const nosh::constant rhs(1.0);

  nosh::function x(mesh);

  nosh::linear_solve(
    matrix, rhs, x,
    {
      {"package", "Belos"},
      {"method", "Pseudo Block CG"},
      {"preconditioner", "MueLu"}
    }
  );

  nosh::write(x, "out.h5m");

  return EXIT_SUCCESS;
}
