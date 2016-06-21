#ifndef NOSH_CONSTANT_HPP
#define NOSH_CONSTANT_HPP

#include <memory>
#include <climits>

#include <Tpetra_Vector.hpp>
#include <Eigen/Dense>

#include "expression.hpp"
#include "mesh.hpp"

namespace nosh {

  class constant: public expression
  {
    public:
      explicit constant(const double val):
        expression(0),
        val_(val)
      {}

      ~constant() override = default;

      double
      operator()(const Eigen::Vector3d & x) const override
      {
        (void) x;
        return val_;
      };

    private:
      const double val_;
  };

} // namespace nosh
#endif // NOSH_CONSTANT_HPP
