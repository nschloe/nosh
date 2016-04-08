#ifndef NOSH_DIRICHLET_BC_CONST_HPP
#define NOSH_DIRICHLET_BC_CONST_HPP

#include "dirichlet_bc.hpp"

namespace nosh {
  class dirichlet_bc_const:
    public dirichlet_bc
  {
    public:
      dirichlet_bc_const(const double val):
        dirichlet_bc({}),
        val_(val)
      {};

      virtual
      ~dirichlet_bc_const()
      {};

      virtual
      bool
      is_inside(const Eigen::Vector3d & x) const
      {
        (void) x;
        return true;
      };

      virtual
      double
      eval(const Eigen::Vector3d & x) const
      {
        (void) x;
        return val_;
      };

    private:
      const double val_;
  };
}
#endif // NOSH_DIRICHLET_BC_CONST_HPP
