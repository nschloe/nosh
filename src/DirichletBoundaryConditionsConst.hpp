#ifndef NOSH_DIRICHLETBOUNDARYCONDITONSCONST_HPP
#define NOSH_DIRICHLETBOUNDARYCONDITONSCONST_HPP

#include "DirichletBoundaryConditions.hpp"

namespace Nosh {
  class DirichletBoundaryConditionsConst:
    public DirichletBoundaryConditions
  {
    public:
      DirichletBoundaryConditionsConst(const double val):
        val_(val)
      {};

      virtual
      ~DirichletBoundaryConditionsConst()
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
