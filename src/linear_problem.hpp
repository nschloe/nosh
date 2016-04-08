#ifndef NOSH_LINEAR_PROBLEM_H
#define NOSH_LINEAR_PROBLEM_H

#include <Teuchos_Tuple.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "helper.hpp"
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
        matrix(mesh->build_graph()),
        rhs(Teuchos::rcp(mesh->map())),
        dbcs_(dbcs),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        fill_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: linear_problem::fill_")),
#endif
        edge_cores_(edge_cores),
        vertex_cores_(vertex_cores),
        boundary_cores_(boundary_cores)
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
#endif
          this->matrix.resumeFill();

          this->matrix.setAllToScalar(0.0);
          this->rhs.putScalar(0.0);

          this->add_edge_contributions_();
          this->add_vertex_contributions_();
          this->add_domain_boundary_contributions_();
          this->apply_dbcs_();

          this->matrix.fillComplete();

          return;
        }
    private:
      void
      add_edge_contributions_()
      {
        for (const auto & edge_core: this->edge_cores_) {
          for (const auto & subdomain_id: edge_core->subdomain_ids) {
            const auto edges = this->mesh_->get_edges(subdomain_id);
            const auto edge_data = this->mesh_->get_edge_data();
            for (const auto edge: edges) {
              const auto verts = this->mesh_->get_vertex_tuple(edge);
              const auto lid = this->mesh_->local_index(edge);
              auto vals = edge_core->eval(
                  this->mesh_->get_coords(verts[0]),
                  this->mesh_->get_coords(verts[1]),
                  edge_data[lid].length,
                  edge_data[lid].covolume
                  );

              const auto & gids = this->mesh_->edge_gids[lid];
              for (int i = 0; i < 2; i++) {
                // Add to matrix
                int num_lhs = this->matrix.sumIntoGlobalValues(
                    gids[i], gids,
                    Teuchos::ArrayView<double>(vals.lhs[i])
                    );
#ifndef NDEBUG
                TEUCHOS_ASSERT_EQUALITY(num_lhs, 2);
#endif
                // Add to rhs
                this->rhs.sumIntoGlobalValue(gids[i], vals.rhs[i]);
              }
            }
          }
        }
      }

      void
      add_vertex_contributions_()
      {
        const auto & control_volumes = this->mesh_->control_volumes();
        const auto c_data = control_volumes->getData();
        for (const auto & vertex_core: this->vertex_cores_) {
          for (const auto & subdomain_id: vertex_core->subdomain_ids) {
            const auto verts = this->mesh_->get_vertices(subdomain_id);
            for (const auto vertex: verts) {
              const auto lid = this->mesh_->local_index(vertex);
              auto val = vertex_core->eval(
                  this->mesh_->get_coords(vertex),
                  c_data[lid]
                  );
              // Add to matrix
              auto gid = this->matrix.getMap()->getGlobalElement(lid);
              int num_lhs = this->matrix.sumIntoGlobalValues(
                  gid,
                  Teuchos::tuple<int>(gid),
                  Teuchos::tuple<double>(val.lhs)
                  );
#ifndef NDEBUG
              TEUCHOS_ASSERT_EQUALITY(num_lhs, 1);
#endif
              // add to rhs
              this->rhs.sumIntoGlobalValue(gid, val.rhs);
            }
          }
        }
      }

      void
      add_domain_boundary_contributions_()
      {
        const auto surfs = this->mesh_->boundary_surface_areas();
        for (const auto boundary_core: this->boundary_cores_) {
          for (const auto & subdomain_id: boundary_core->subdomain_ids) {
            const auto verts = this->mesh_->get_vertices(subdomain_id);
            for (const auto vert: verts) {
              const int lid = this->mesh_->local_index(vert);
              auto val = boundary_core->eval(
                  this->mesh_->get_coords(vert),
                  surfs[lid]
                  );
              // Add to matrix
              const auto gid = this->matrix.getMap()->getGlobalElement(lid);
              int num_lhs = this->matrix.sumIntoGlobalValues(
                  gid,
                  Teuchos::tuple<int>(gid),
                  Teuchos::tuple<double>(val.lhs)
                  );
#ifndef NDEBUG
              TEUCHOS_ASSERT_EQUALITY(num_lhs, 1);
#endif
              // add to rhs
              this->rhs.sumIntoGlobalValue(gid, val.rhs);
            }
          }
        }
      }

      // apply dirichlet boundary conditions
      void
      apply_dbcs_()
      {
        for (const auto & bc: this->dbcs_) {
          for (const auto & subdomain_id: bc->subdomain_ids) {
            const auto verts = this->mesh_->get_vertices(subdomain_id);
            for (const auto & vertex: verts) {
              // eliminate the row in A
              const auto gid = this->mesh_->gid(vertex);
              size_t num = this->matrix.getNumEntriesInGlobalRow(gid);
              // It shouldn't actually happen that the specified global row
              // does not belong to this graph.
              // TODO find out if why we need this, fix the underlying issue,
              // make this a TEUCHOS_TEST_*
              if (num != Teuchos::OrdinalTraits<size_t>::invalid()) {
                std::vector<int> cols(num);
                std::vector<double> vals(num);
                this->matrix.getGlobalRowCopy(gid, cols, vals, num);
                // set vals to 0
                std::fill(vals.begin(), vals.end(), 0.0);
                // set diagonal entry to 1
                auto it = std::find(cols.begin(), cols.end(), gid);
                TEUCHOS_TEST_FOR_EXCEPT_MSG(
                    it == cols.end(),
                    "Matrix has no main diagonal entry."
                    );
                int pos = it - cols.begin();
                // set diagonal entry to 1
                vals[pos] = 1.0;
                this->matrix.replaceGlobalValues(gid, cols, vals);
                // set rhs
                const auto coord = this->mesh_->get_coords(vertex);
                this->rhs.replaceGlobalValue(gid, bc->eval(coord));
              }
            }
          }
        }
      }

    private:
      const std::shared_ptr<const nosh::mesh> mesh_;

    public:
      Tpetra::CrsMatrix<double,int,int> matrix;
      Tpetra::Vector<double,int,int> rhs;

    private:
      // TODO const?
      const std::set<std::shared_ptr<const nosh::dirichlet_bc>> dbcs_;
#ifdef NOSH_TEUCHOS_TIME_MONITOR
      const Teuchos::RCP<Teuchos::Time> fill_time_;
#endif
      const std::set<std::shared_ptr<const edge_core>> edge_cores_;
      const std::set<std::shared_ptr<const vertex_core>> vertex_cores_;
      const std::set<std::shared_ptr<const boundary_core>> boundary_cores_;
  };
} // namespace nosh

#endif // NOSH_LINEAR_PROBLEM_H
