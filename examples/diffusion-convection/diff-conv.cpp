#include "diff-conv.hpp"
#include <nosh.hpp>
#include <memory>

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  diff_conv::laplace matrix(mesh);

  const diff_conv::f rhs;

  nosh::function x(mesh);

  nosh::scaled_linear_solve(
      matrix, rhs, x,
      {
        {"package", "Belos"},
        {"method", "Pseudo Block GMRES"}
      });

  nosh::write(x, "out.h5m");

  return EXIT_SUCCESS;
}
