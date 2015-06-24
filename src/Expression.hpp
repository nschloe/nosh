#ifndef NOSH_EXPRESSION_HPP
#define NOSH_EXPRESSION_HPP

#include <memory>
#include <climits>

#include <Tpetra_Vector.hpp>
#include <Eigen/Dense>

#include "Mesh.hpp"

namespace Nosh {

  class Expression
  {
    public:
      Expression(const unsigned int _degree = UINT_MAX):
        degree(_degree)
      {}

      virtual
      ~Expression()
      {}

      virtual
      double
      eval(const Eigen::Vector3d & x) const = 0;

    public:
      const unsigned int degree;
  };

  // helper functions
  std::shared_ptr<Tpetra::Vector<double,int,int>>
  integrateOverControlVolumes(
      const Nosh::Expression & expr,
      const Nosh::Mesh & mesh
      );

}
#endif // NOSH_EXPRESSION_HPP
