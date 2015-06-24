#include "Poisson.hpp"
#include <nosh.hpp>
#include <memory>

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  //const auto mesh = Nosh::read("rectangle.e");
  const auto mesh = Nosh::read("pacman.e");
  //const auto mesh = Nosh::read("cubesmall.e");
  //const auto mesh = Nosh::read("brick-w-hole.e");

  const auto bc1 = std::make_shared<Poisson::BC1>();
  const auto bc2 = std::make_shared<Poisson::BC2>();
  const auto bc3 = std::make_shared<Poisson::BC3>();

  //const Poisson::A A(mesh, {bc1, bc2});
  const Poisson::A A(mesh, {bc3});

  //Nosh::Constant f(1.0);
  Poisson::F f;

  Nosh::Function x(mesh);
  x.putScalar(0.0);

  Nosh::linearSolve(A, f, x);

  Nosh::write(x, "out.e");

  return EXIT_SUCCESS;
}