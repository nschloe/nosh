#ifndef NOSH_CONSTANT_HPP
#define NOSH_CONSTANT_HPP

#include <memory>
#include <climits>

#include <Tpetra_Vector.hpp>
#include <Eigen/Dense>

#include "mesh.hpp"

namespace nosh {

  class constant: public expression
  {
    public:
      constant(const double val):
        expression(0),
        val_(val)
      {}

      virtual
      ~constant()
      {}

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
#endif // NOSH_CONSTANT_HPP
