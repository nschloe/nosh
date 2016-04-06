#ifndef NOSH_MATRIX_HPP
#define NOSH_MATRIX_HPP

#include <memory>
#include <set>

#include <Tpetra_CrsMatrix.hpp>

#include "dirichlet_bc.hpp"
#include "mesh.hpp"

namespace nosh {
  class matrix:
    public Tpetra::CrsMatrix<double,int,int>
  {
    public:
      matrix(
          const std::shared_ptr<const nosh::mesh> & _mesh,
          const std::set<std::shared_ptr<const nosh::dirichlet_bc>> & _bcs
          ):
        Tpetra::CrsMatrix<double,int,int>(_mesh->build_graph()),
        mesh(_mesh),
        bcs(_bcs)
    {
    };

    virtual
    ~matrix()
    {};

    public:
    const std::shared_ptr<const nosh::mesh> mesh;
    const std::set<std::shared_ptr<const nosh::dirichlet_bc>> bcs;

    protected:

    void
    apply_bcs_();
  };
}
#endif // NOSH_MATRIX_HPP
