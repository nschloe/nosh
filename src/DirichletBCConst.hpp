#ifndef NOSH_DIRICHLETBOUNDARYCONDITONSCONST_HPP
#define NOSH_DIRICHLETBOUNDARYCONDITONSCONST_HPP

#include "DirichletBC.hpp"

namespace Nosh {
  class DirichletBCConst:
    public DirichletBC
  {
    public:
      DirichletBCConst(const double val):
        val_(val)
      {};

      virtual
      ~DirichletBCConst()
      {};

      virtual
      bool
      isInside(const Eigen::Vector3d & x) const
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
#endif // NOSH_DIRICHLETBOUNDARYCONDITONSCONST_HPP
