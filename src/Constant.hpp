#ifndef NOSH_CONSTANT_HPP
#define NOSH_CONSTANT_HPP

#include <memory>
#include <climits>

#include <Tpetra_Vector.hpp>
#include <Eigen/Dense>

#include "Mesh.hpp"

namespace Nosh {

  class Constant: public Expression
  {
    public:
      Constant(const double val):
        Expression(0),
        val_(val)
      {}

      virtual
      ~Constant()
      {}

      virtual
      double
      eval(const Eigen::Vector3d & x) const
      {
        (void) x;
        return value_;
      };

    private:
      const double val_;
  };

}
#endif // NOSH_CONSTANT_HPP
