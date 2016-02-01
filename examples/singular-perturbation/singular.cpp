#include "singular.hpp"
#include <nosh.hpp>
#include <memory>

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  const auto bc1 = std::make_shared<singular::bc1>();

  singular::singular matrix(mesh, {bc1});

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
