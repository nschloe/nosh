#ifndef NOSH_LINEAR_PROBLEM_H
#define NOSH_LINEAR_PROBLEM_H

#include <Teuchos_Tuple.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "helper.hpp"
#include "fvm_helpers.hpp"
#include "matrix.hpp"
#include "edge_core.hpp"
#include "vertex_core.hpp"
#include "boundary_core.hpp"
#include "mesh.hpp"

namespace nosh
{
  class linear_problem
  {
    public:
      linear_problem(
          const std::shared_ptr<const nosh::mesh> & mesh,
          const std::set<std::shared_ptr<const edge_core>> & edge_cores,
          const std::set<std::shared_ptr<const vertex_core>> & vertex_cores,
          const std::set<std::shared_ptr<const boundary_core>> & boundary_cores,
          const std::set<std::shared_ptr<const dirichlet_bc>> & dbcs
          ) :
        mesh_(mesh),
        matrix(std::make_shared<Tpetra::CrsMatrix<double,int,int>>(mesh->build_graph())),
        rhs(std::make_shared<Tpetra::Vector<double,int,int>>(Teuchos::rcp(mesh->map()))),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        fill_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: linear_problem::fill_")),
#endif
        edge_cores_(edge_cores),
        vertex_cores_(vertex_cores),
        boundary_cores_(boundary_cores),
        dbcs_(dbcs)
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
          TEUCHOS_ASSERT(this->mesh_);
          TEUCHOS_ASSERT(this->matrix);
          TEUCHOS_ASSERT(this->rhs);
#endif
          this->matrix->resumeFill();

          this->matrix->setAllToScalar(0.0);
          this->rhs->putScalar(0.0);

          fvm_helpers::add_edge_contributions(
              this->mesh_,
              this->edge_cores_,
              this->matrix,
              this->rhs
              );
          fvm_helpers::add_vertex_contributions(
              this->mesh_,
              this->vertex_cores_,
              this->matrix,
              this->rhs
              );
          fvm_helpers::add_domain_boundary_contributions(
              this->mesh_,
              this->boundary_cores_,
              this->matrix,
              this->rhs
              );
          fvm_helpers::apply_dbcs(
              this->mesh_,
              this->dbcs_,
              this->matrix,
              this->rhs
              );

          this->matrix->fillComplete();

          return;
        }

    private:
      const std::shared_ptr<const nosh::mesh> mesh_;

    public:
      const std::shared_ptr<Tpetra::CrsMatrix<double,int,int>> matrix;
      const std::shared_ptr<Tpetra::Vector<double,int,int>> rhs;

    private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
      const Teuchos::RCP<Teuchos::Time> fill_time_;
#endif
      const std::set<std::shared_ptr<const edge_core>> edge_cores_;
      const std::set<std::shared_ptr<const vertex_core>> vertex_cores_;
      const std::set<std::shared_ptr<const boundary_core>> boundary_cores_;
      const std::set<std::shared_ptr<const nosh::dirichlet_bc>> dbcs_;
  };
} // namespace nosh

#endif // NOSH_LINEAR_PROBLEM_H
