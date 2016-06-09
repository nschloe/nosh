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
      {"package", std::string("Belos")},
      {"method", std::string("Pseudo Block CG")},
      {"preconditioner", std::string("MueLu")}
    }
  );

  nosh::write(x, "out.h5m");

  return EXIT_SUCCESS;
}
