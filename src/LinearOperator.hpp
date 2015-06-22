#ifndef NOSH_LINEAROPERATOR_HPP
#define NOSH_LINEAROPERATOR_HPP

#include <memory>

#include <Tpetra_CrsMatrix.hpp>

#include "DirichletBoundaryConditions.hpp"
#include "Mesh.hpp"

namespace Nosh {
  class LinearOperator:
    public Tpetra::CrsMatrix<double,int,int>
  {
    public:
      LinearOperator(
          const std::shared_ptr<const Nosh::Mesh> & _mesh,
          const std::shared_ptr<const Nosh::DirichletBoundaryConditions> & _bcs
          ):
        Tpetra::CrsMatrix<double,int,int>(_mesh->buildGraph()),
        mesh(_mesh),
        bcs(_bcs)
    {
    };

    virtual
    ~LinearOperator()
    {};

    public:
    const std::shared_ptr<const Nosh::Mesh> mesh;
    const std::shared_ptr<const Nosh::DirichletBoundaryConditions> bcs;
  };
}
#endif // NOSH_LINEAROPERATOR_HPP
