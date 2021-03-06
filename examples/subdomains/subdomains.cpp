#include "subdomains.hpp"
#include <nosh.hpp>
#include <memory>

using dict = std::map<std::string, boost::any>;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  const auto mesh = nosh::read("pacman.h5m");

  mesh->mark_subdomains({
      std::make_shared<subdomains::d1>(),
      std::make_shared<subdomains::d2>()
      });

  subdomains::p problem(mesh);

  nosh::function x(mesh);
  x.putScalar(0.0);

  mikado::linear_solve(
      *(problem.matrix), *(problem.rhs), x,
      {
        {"package", std::string("Belos")},
        //,{"method", std::string("Pseudo Block GMRES")},
        {"method", std::string("Pseudo Block CG")},
        {"parameters", dict{
          {"Convergence Tolerance", 1.0e-10},
          {"Output Frequency", 1},
          {"Output Style", 1},
          {"Verbosity", 33}
        }},
        {"preconditioner", std::string("MueLu")}
        // {"preconditioner matrix", M},
        // {"preconditioner parameters", dict{
        //   {"cycle type", std::string("V")}
        // }}
      }
      );

  nosh::write(x, "out.h5m");

  return EXIT_SUCCESS;
}
