#include <nosh.hpp>

#include <Teuchos_StandardCatchMacros.hpp>

namespace different_materials {
  class eps: public nosh::expression
  {
    public:
      eps(): nosh::expression(0) {}
      virtual ~eps() {}

      virtual double operator()(const Eigen::Vector3d & x) const
      {
        if (x[0] > 0.0) {
          return 3.0;
        } else {
          return 1.0;
        }
      };
  };
}

#include "different-materials.hpp"

using list = std::map<std::string, boost::any>;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  auto out = Teuchos::VerboseObjectBase::getDefaultOStream();

  bool success = true;
  try {
  const auto mesh = nosh::read("pacman.h5m");

  different_materials::problem p(mesh);

  nosh::function x(mesh);
  x.putScalar(0.0);

  //nosh::scaled_linear_solve(
  nosh::linear_solve(
      p, x,
      {
        {"package", "Belos"},
        {"method", "Pseudo Block GMRES"},
        {"parameters", list{
          {"Convergence Tolerance", 1.0e-10},
          {"Output Frequency", 1},
          {"Output Style", 1},
          {"Verbosity", 33}
        }},
        {"preconditioner", "MueLu"}
      }
      );

  nosh::write(x, "out.h5m");
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

  return EXIT_SUCCESS;
}
