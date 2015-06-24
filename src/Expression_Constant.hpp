#ifndef NOSH_CONSTANT_HPP
#define NOSH_CONSTANT_HPP

#include "Expression.hpp"

namespace Nosh {

  class Constant:
    public Expression
  {
    public:
      Constant(const double alpha):
        Expression(0),
        alpha_(alpha)
      {}

      virtual
      ~Constant()
      {}

      virtual
      double
      eval(const Eigen::Vector3d & x) const
      {
        (void) x;
        return alpha_;
      };

    private:
      const double alpha_;
  };

}
#endif // NOSH_CONSTANT_HPP
