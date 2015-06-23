#ifndef NOSH_DIRICHLETBOUNDARYCONDITONS_HPP
#define NOSH_DIRICHLETBOUNDARYCONDITONS_HPP

#include <Eigen/Dense>

namespace Nosh {
  class DirichletBC
  {
    public:
      DirichletBC()
      {};

      virtual
      ~DirichletBC()
      {};

      virtual
      bool
      isInside(const Eigen::Vector3d & x) const = 0;

      virtual
      double
      eval(const Eigen::Vector3d & x) const = 0;
  };
}
#endif // NOSH_DIRICHLETBOUNDARYCONDITONS_HPP
