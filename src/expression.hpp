#ifndef NOSH_EXPRESSION_HPP
#define NOSH_EXPRESSION_HPP

#include <memory>
#include <climits>

#include <Tpetra_Vector.hpp>
#include <Eigen/Dense>

#include "mesh.hpp"

namespace nosh {

  class expression
  {
    public:
      expression(const unsigned int _degree = UINT_MAX):
        degree(_degree)
      {}

      virtual
      ~expression()
      {}

      virtual
      double
      eval(const Eigen::Vector3d & x) const = 0;

    public:
      const unsigned int degree;
  };

  // helper functions
  std::shared_ptr<Tpetra::Vector<double,int,int>>
  integrate_over_control_volumes(
      const nosh::expression & expr,
      const nosh::mesh & mesh
      );

}
#endif // NOSH_EXPRESSION_HPP
