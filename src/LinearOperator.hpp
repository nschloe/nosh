#ifndef NOSH_LINEAROPERATOR_HPP
#define NOSH_LINEAROPERATOR_HPP

#include <memory>
#include <set>

#include <Tpetra_CrsMatrix.hpp>

#include "DirichletBC.hpp"
#include "Mesh.hpp"

namespace Nosh {
  class LinearOperator:
    public Tpetra::CrsMatrix<double,int,int>
  {
    public:
      LinearOperator(
          const std::shared_ptr<const Nosh::Mesh> & _mesh,
          const std::set<std::shared_ptr<const Nosh::DirichletBC>> & _bcs
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
    const std::set<std::shared_ptr<const Nosh::DirichletBC>> bcs;

    protected:

    void
    applyBcs_();
  };
}
#endif // NOSH_LINEAROPERATOR_HPP
