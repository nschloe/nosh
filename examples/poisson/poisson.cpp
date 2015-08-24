#include "poisson.hpp"
#include <nosh.hpp>
#include <memory>

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.e");

  const auto bc1 = std::make_shared<poisson::bc1>();
  const auto bc2 = std::make_shared<poisson::bc2>();

  poisson::laplace matrix(mesh, {bc1, bc2});

  const poisson::f rhs;

  nosh::function x(mesh);

  nosh::scaled_linear_solve(
      matrix, rhs, x
      );

  nosh::write(x, "out.e");

  return EXIT_SUCCESS;
}
