#include "singular.hpp"
#include <nosh.hpp>

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  const singular::singular problem(mesh);

  nosh::function x(mesh);

  nosh::linear_solve(
    problem, x,
    {
      {"package", "Belos"},
      {"method", "Pseudo Block CG"},
      {"preconditioner", "MueLu"}
    }
  );

  nosh::write(x, "out.h5m");

  return EXIT_SUCCESS;
}
