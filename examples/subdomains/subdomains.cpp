#include "subdomains.hpp"
#include <nosh.hpp>
#include <memory>

using list = std::map<std::string, boost::any>;
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

  //nosh::scaled_linear_solve(
  nosh::linear_solve(
      problem, x,
      {
        {"package", "Belos"},
        //,{"method", "Pseudo Block GMRES"},
        {"method", "Pseudo Block CG"},
        {"parameters", list{
          {"Convergence Tolerance", 1.0e-10},
          {"Output Frequency", 1},
          {"Output Style", 1},
          {"Verbosity", 33}
        }},
        {"preconditioner", "MueLu"}
        // {"preconditioner matrix", M},
        // {"preconditioner parameters", list{
        //   {"cycle type", "V"}
        // }}
      }
      );

  nosh::write(x, "out.h5m");

  return EXIT_SUCCESS;
}
