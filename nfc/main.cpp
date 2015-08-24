#include "poisson.hpp"
#include <nosh.hpp>
#include <memory>

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  // auto p = Teuchos::rcp(new Teuchos::ParameterList());
  // nosh::stdmap2teuchoslist(nosh::defaultLinearSolverParams(), *p);
  // std::cout << *p << std::endl;

  //const auto mesh = nosh::read("rectangle.e");
  const auto mesh = nosh::read("pacman.e");
  //const auto mesh = nosh::read("cubesmall.e");
  //const auto mesh = nosh::read("brick-w-hole.e");

  const auto bc1 = std::make_shared<poisson::bc1>();
  const auto bc2 = std::make_shared<poisson::bc2>();
  const auto bc3 = std::make_shared<poisson::bc3>();

  //poisson::A A(mesh, {bc1, bc2});
  //poisson::A2 matrix(mesh, {bc3});
  poisson::i matrix(mesh, {bc3});

  nosh::constant rhs(1.0);
  //const poisson::f rhs();

  nosh::function x(mesh);
  x.putScalar(0.0);

  nosh::scaled_linear_solve(
      matrix, rhs, x
      // {
      // {"method", "GMRES"},
      // {"parameters", {
      //   // ...
      // }}
      // }
      );

  nosh::write(x, "out.e");

  return EXIT_SUCCESS;
}
