#include "diff-conv.hpp"

#include <nosh.hpp>
#include <mikado.hpp>
#include <memory>

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  diff_conv::dc problem(mesh);

  nosh::function x(mesh);

  // TODO
  // nosh::scaled_linear_solve(
  mikado::linear_solve(
      problem, x,
      {
        {"package", "Belos"},
        {"method", "Pseudo Block GMRES"}
      });

  nosh::write(x, "out.h5m");

  return EXIT_SUCCESS;
}
