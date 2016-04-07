#ifndef NOSH_FVMMATRIX_H
#define NOSH_FVMMATRIX_H

#include <Teuchos_Tuple.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "helper.hpp"
#include "matrix.hpp"
#include "matrix_core.hpp"
#include "mesh.hpp"

namespace nosh
{
  class linear_problem
  {
    public:
      linear_problem(
          const std::shared_ptr<const nosh::mesh> & mesh,
          const std::set<std::shared_ptr<const matrix_core>> & matrix_cores,
          const std::set<std::shared_ptr<const nosh::dirichlet_bc>> & dbcs
          ) :
        mesh_(mesh),
        matrix(mesh->build_graph()),
        rhs(Teuchos::rcp(mesh->map())),
        dbcs_(dbcs),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        fill_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: linear_problem::fill_"))
#endif
        // ,matrix_cores_(matrix_cores), TODO
        ,edge_core_associations_(this->associate_edge_cores(matrix_cores))
        ,vertex_core_associations_(this->associate_vertex_cores(matrix_cores))
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
      // Create a vector core_idx that associates exactly one core with each
      // edge. This core will be the one used for building the matrix at that
      // edge.
      std::vector<std::shared_ptr<const matrix_core>>
        associate_edge_cores(
          const std::set<std::shared_ptr<const matrix_core>> & matrix_cores
          )
      {
        const std::vector<edge> edges = this->mesh_->my_edges();
        std::vector<std::shared_ptr<const matrix_core>> core_idx(edges.size());
        for (size_t k = 0; k < edges.size(); k++) {
          const auto x0 = this->mesh_->get_coords(std::get<0>(edges[k]));
          const auto x1 = this->mesh_->get_coords(std::get<1>(edges[k]));
          const auto mp = 0.5 * (x0 + x1);
          bool is_assigned = false;
          for (const auto mc: matrix_cores) {
            if (mc->is_inside(mp)) {
              is_assigned = true;
              core_idx[k] = mc;
              break;
            }
          }
          TEUCHOS_ASSERT(is_assigned);
        }
        return core_idx;
      }

      std::vector<std::shared_ptr<const matrix_core>>
        associate_vertex_cores(
          const std::set<std::shared_ptr<const matrix_core>> & matrix_cores
          )
      {
        const auto & owned_nodes = this->mesh_->get_owned_nodes();
        std::vector<std::shared_ptr<const matrix_core>> core_idx(owned_nodes.size());
        for (size_t k = 0; k < owned_nodes.size(); k++) {
          const auto x = this->mesh_->get_coords(owned_nodes[k]);
          bool is_assigned = false;
          for (const auto mc: matrix_cores) {
            if (mc->is_inside(x)) {
              is_assigned = true;
              core_idx[k] = mc;
              break;
            }
          }
          TEUCHOS_ASSERT(is_assigned);
        }
        return core_idx;
      }

      void
      add_edge_contributions_()
      {
        const std::vector<edge> edges = this->mesh_->my_edges();
        const auto edge_data = this->mesh_->get_edge_data();
        for (size_t k = 0; k < edges.size(); k++) {
          auto vals = this->edge_core_associations_[k]->edge_contrib(
              this->mesh_->get_coords(std::get<0>(edges[k])),
              this->mesh_->get_coords(std::get<1>(edges[k])),
              edge_data[k].length,
              edge_data[k].covolume
              );

          const Teuchos::Tuple<int,2> & gids = this->mesh_->edge_gids[k];
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

      void
      add_vertex_contributions_()
      {
        const auto & control_volumes = this->mesh_->control_volumes();
        const auto & owned_nodes = this->mesh_->get_owned_nodes();
        auto c_data = control_volumes->getData();
#ifndef NDEBUG
        TEUCHOS_ASSERT_EQUALITY(c_data.size(), owned_nodes.size());
#endif
        for (int k = 0; k < c_data.size(); k++) {
          auto val = this->vertex_core_associations_[k]->vertex_contrib(
              this->mesh_->get_coords(owned_nodes[k]),
              c_data[k]
              );
          // Add to matrix
          auto gid = this->matrix.getMap()->getGlobalElement(k);
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

      void
      add_domain_boundary_contributions_()
      {
        const auto verts = this->mesh_->boundary_vertices();
        const auto surfs = this->mesh_->boundary_surface_areas();
        const auto & owned_nodes = this->mesh_->get_owned_nodes();
#ifndef NDEBUG
        TEUCHOS_ASSERT_EQUALITY(verts.size(), surfs.size());
#endif
        for (size_t k = 0; k < verts.size(); k++) {
          const int lid = this->mesh_->local_index(verts[k]);
          auto val = this->vertex_core_associations_[lid]->domain_boundary_contrib(
              this->mesh_->get_coords(owned_nodes[lid]),
              surfs[k]
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


      // apply dirichlet boundary conditions
      void
      apply_dbcs_()
      {
        const auto boundary_vertices = this->mesh_->boundary_vertices();
        for (const auto boundary_vertex: boundary_vertices) {
          // check if any of the boundary conditions kicks in
          const auto coord = this->mesh_->get_coords(boundary_vertex);
          int count = 0;
          std::shared_ptr<const nosh::dirichlet_bc> active_bc;
          for (const auto & bc: this->dbcs_) {
            if (bc->is_inside(coord)) {
              active_bc = bc;
              count++;
            }
          }
          TEUCHOS_TEST_FOR_EXCEPT_MSG(
              count > 1,
              "More than one active boundary conditions. Abort."
              );
          if (count == 1) {
            // eliminate the row in A
            const auto gid = this->mesh_->gid(boundary_vertex);
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
              this->rhs.replaceGlobalValue(gid, active_bc->eval(coord));
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
      std::vector<std::shared_ptr<const matrix_core>> edge_core_associations_;
      std::vector<std::shared_ptr<const matrix_core>> vertex_core_associations_;
  };
} // namespace nosh

#endif // NOSH_FVMMATRIX_H
