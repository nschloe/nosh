#ifndef NOSH_DIRICHLET_BC_HPP
#define NOSH_DIRICHLET_BC_HPP

#include <Eigen/Dense>

namespace nosh {
  class dirichlet_bc
  {
    public:
      dirichlet_bc()
      {};

      virtual
      ~dirichlet_bc()
      {};

      virtual
      bool
      is_inside(const Eigen::Vector3d & x) const = 0;

      virtual
      double
      eval(const Eigen::Vector3d & x) const = 0;
  };
}
#endif // NOSH_DIRICHLET_BC_HPP
