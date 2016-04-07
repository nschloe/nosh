#ifndef NOSH_DIRICHLET_BC_HPP
#define NOSH_DIRICHLET_BC_HPP

#include <Eigen/Dense>

namespace nosh {
  class dirichlet_bc
  {
    public:
      dirichlet_bc(const std::set<std::string> & _subdomain_ids):
        subdomain_ids(_subdomain_ids)
      {};

      virtual
      ~dirichlet_bc()
      {};

      virtual
      double
      eval(const Eigen::Vector3d & x) const = 0;

    public:
      const std::set<std::string> subdomain_ids;
  };
}
#endif // NOSH_DIRICHLET_BC_HPP
