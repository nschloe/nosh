#include "Poisson.hpp"
#include <nosh.hpp>
#include <memory>

class BC1: public Nosh::DirichletBC {
  virtual bool
  isInside(const Eigen::Vector3d & x) const {
    return x[0] < 0.0;
  }
  virtual double
  eval(const Eigen::Vector3d & x) const {
    (void) x;
    return 0.0;
  }
};

class BC2: public Nosh::DirichletBC {
  virtual bool
  isInside(const Eigen::Vector3d & x) const {
    return x[0] >= 0.0;
  }
  virtual double
  eval(const Eigen::Vector3d & x) const {
    (void) x;
    return 1.0;
  }
};

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  //const auto mesh = Nosh::read("rectangle.e");
  const auto mesh = Nosh::read("pacman.e");
  //const auto mesh = Nosh::read("cubesmall.e");
  //const auto mesh = Nosh::read("brick-w-hole.e");

  const auto bc1 = std::make_shared<BC1>();
  const auto bc2 = std::make_shared<BC2>();

  const Poisson::A A(mesh, {bc1, bc2});

  Nosh::Function f(mesh);
  f.putScalar(0.0);

  Nosh::Function x(mesh);
  x.putScalar(0.0);

  Nosh::linearSolve(A, f, x);

  Nosh::write(x, "out.e");

  return EXIT_SUCCESS;
}
