#ifndef NOSH_LINEAR_PROBLEM_H
#define NOSH_LINEAR_PROBLEM_H

#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "fvm_matrix.hpp"
#include "matrix_core_edge.hpp"
#include "matrix_core_vertex.hpp"
#include "matrix_core_boundary.hpp"
#include "mesh.hpp"

namespace nosh
{
  class linear_problem
  {
    public:
      linear_problem(
          const std::shared_ptr<const nosh::mesh> & mesh,
          const std::set<std::shared_ptr<const matrix_core_edge>> & matrix_core_edges,
          const std::set<std::shared_ptr<const matrix_core_vertex>> & matrix_core_vertexs,
          const std::set<std::shared_ptr<const matrix_core_boundary>> & matrix_core_boundarys,
          const std::set<std::shared_ptr<const matrix_core_dirichlet>> & dbcs
          ) :
        mesh_(mesh),
        matrix(std::make_shared<nosh::fvm_matrix>(mesh, matrix_core_edges, matrix_core_vertexs, matrix_core_boundarys, dbcs)),
        rhs(std::make_shared<Tpetra::Vector<double,int,int>>(Teuchos::rcp(mesh->map())))
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        ,fill_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: linear_problem::fill_"))
#endif
        {
        }

      virtual
      ~linear_problem()
      {};

    protected:
      void
        fill_()
        {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor tm(*fill_time_);
#endif
#ifndef NDEBUG
          TEUCHOS_ASSERT(this->matrix);
          TEUCHOS_ASSERT(this->rhs);
#endif
          // Fill matrix and rhs
          this->matrix->fill(rhs);

          return;
        }

    private:
      const std::shared_ptr<const nosh::mesh> mesh_;

    public:
      const std::shared_ptr<nosh::fvm_matrix> matrix;
      const std::shared_ptr<Tpetra::Vector<double,int,int>> rhs;

    private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
      const Teuchos::RCP<Teuchos::Time> fill_time_;
#endif
  };
} // namespace nosh

#endif // NOSH_LINEAR_PROBLEM_H
