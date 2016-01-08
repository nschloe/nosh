#include "poisson.hpp"
#include <nosh.hpp>
#include <memory>

#include <Teuchos_StandardCatchMacros.hpp>

using list = std::map<std::string, boost::any>;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  auto out = Teuchos::VerboseObjectBase::getDefaultOStream();

  bool success = true;
  try {

  //const auto mesh = nosh::read("r2.h5m");
  const auto mesh = nosh::read("pacman2.h5m");

  const auto bc1 = std::make_shared<poisson::bc1>();
  const auto bc2 = std::make_shared<poisson::bc2>();

  poisson::laplace matrix(mesh, {bc1, bc2});

  const poisson::f rhs;

  nosh::function x(mesh);
  x.putScalar(0.0);

  nosh::scaled_linear_solve(
      matrix, rhs, x,
      {
#if 0
        {"package", "MueLu"}
#endif
#if 1
        // Check
        // https://trilinos.org/docs/dev/packages/amesos2/doc/html/group__amesos2__solver__parameters.html
        // for more options.
        {"package", "Amesos2"},
        {"parameters", list{
          {"Trans", "NOTRANS"},
          {"ColPerm", "COLAMD"}
        }}
#endif
#if 0
        {"package", "Belos"}
        ,{"method", "Pseudo Block GMRES"}
        ,{"parameters", list{
          {"Convergence Tolerance", 1.0e-10},
          {"Output Frequency", 1},
          {"Output Style", 1},
          {"Verbosity", 33}
        }}
        ,{"preconditioner", "MueLu"},
        {"preconditioner parameters", list{
        }}
#endif
      }
      );

  nosh::write(x, "out.h5m");
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success);

  return EXIT_SUCCESS;
}
