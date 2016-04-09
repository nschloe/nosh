#ifndef NOSH_FVM_MATRIX_H
#define NOSH_FVM_MATRIX_H

#include <Tpetra_CrsMatrix.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "edge_core.hpp"
#include "vertex_core.hpp"
#include "boundary_core.hpp"
#include "dirichlet_bc.hpp"
#include "mesh.hpp"

namespace nosh
{
  class fvm_matrix:
    public Tpetra::CrsMatrix<double,int,int>
  {
    public:
      fvm_matrix(
          const std::shared_ptr<const nosh::mesh> & _mesh,
          const std::set<std::shared_ptr<const edge_core>> & edge_cores,
          const std::set<std::shared_ptr<const vertex_core>> & vertex_cores,
          const std::set<std::shared_ptr<const boundary_core>> & boundary_cores,
          const std::set<std::shared_ptr<const dirichlet_bc>> & dbcs
          ) :
        Tpetra::CrsMatrix<double,int,int>(_mesh->build_graph()),
        mesh(_mesh),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        fill_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: fvm_matrix::fill_")),
#endif
        edge_cores_(edge_cores),
        vertex_cores_(vertex_cores),
        boundary_cores_(boundary_cores),
        dbcs_(dbcs)
        {
        }

      virtual
      ~fvm_matrix()
      {};

    public:
      void
        fill(
          const std::shared_ptr<Tpetra::Vector<double,int,int>> & rhs = nullptr
          )
        {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor tm(*fill_time_);
#endif
#ifndef NDEBUG
          TEUCHOS_ASSERT(this->mesh);
#endif
          this->resumeFill();

          this->setAllToScalar(0.0);

          this->add_edge_contributions_(rhs);
          this->add_vertex_contributions_(rhs);
          this->add_domain_boundary_contributions_(rhs);
          this->apply_dbcs_(rhs);

          this->fillComplete();

          return;
        }

    private:

      void
      add_edge_contributions_(
          const std::shared_ptr<Tpetra::Vector<double,int,int>> & rhs
          )
      {
        for (const auto & edge_core: this->edge_cores_) {
          for (const auto & subdomain_id: edge_core->subdomain_ids) {
            const auto edge_data = this->mesh->get_edge_data();

            // this->meshset interior edges
            const auto edges = this->mesh->get_edges(subdomain_id);
            for (const auto edge: edges) {
              const auto verts = this->mesh->get_vertex_tuple(edge);
              const auto lid = this->mesh->local_index(edge);
              auto vals = edge_core->eval(
                  this->mesh->get_coords(verts[0]),
                  this->mesh->get_coords(verts[1]),
                  edge_data[lid].length,
                  edge_data[lid].covolume
                  );

              const auto & gids = this->mesh->edge_gids[lid];
              for (int i = 0; i < 2; i++) {
                // Add to matrix
                const int num_lhs = this->sumIntoGlobalValues(
                    gids[i], gids,
                    Teuchos::ArrayView<double>(vals.lhs[i])
                    );
#ifndef NDEBUG
                TEUCHOS_ASSERT_EQUALITY(num_lhs, 2);
#endif
                if (rhs) {
                  // Add to rhs
                  rhs->sumIntoGlobalValue(gids[i], vals.rhs[i]);
                }
              }
            }

            // this->meshset boundary edges
            const auto half_edges = this->mesh->get_edges(
                subdomain_id + "_halfedges"
                );
            for (const auto edge: half_edges) {
              const auto verts = this->mesh->get_vertex_tuple(edge);

              // check which one of the two verts is in this->meshset
              int i;
              if (this->mesh->contains(subdomain_id, {verts[0]})) {
                i = 0;
              } else if (this->mesh->contains(subdomain_id, {verts[1]})) {
                i = 1;
              } else {
                TEUCHOS_TEST_FOR_EXCEPT_MSG(
                    true,
                    "Neither of the two edge vertices is contained in the subdomain."
                    );
              }

              const auto lid = this->mesh->local_index(edge);
              auto vals = edge_core->eval(
                  this->mesh->get_coords(verts[0]),
                  this->mesh->get_coords(verts[1]),
                  edge_data[lid].length,
                  edge_data[lid].covolume
                  );

              const auto & gids = this->mesh->edge_gids[lid];
              // Add to matrix
              int num_lhs = this->sumIntoGlobalValues(
                  gids[i], gids,
                  Teuchos::ArrayView<double>(vals.lhs[i])
                  );
#ifndef NDEBUG
              TEUCHOS_ASSERT_EQUALITY(num_lhs, 2);
#endif
              if (rhs) {
                // Add to rhs
                rhs->sumIntoGlobalValue(gids[i], vals.rhs[i]);
              }
            }
          }
        }
      }

      void
      add_vertex_contributions_(
          const std::shared_ptr<Tpetra::Vector<double,int,int>> & rhs
          )
      {
        const auto & control_volumes = this->mesh->control_volumes();
        const auto c_data = control_volumes->getData();
        for (const auto & vertex_core: this->vertex_cores_) {
          for (const auto & subdomain_id: vertex_core->subdomain_ids) {
            const auto verts = this->mesh->get_vertices(subdomain_id);
            for (const auto & vertex: verts) {
              const auto lid = this->mesh->local_index(vertex);
              const auto val = vertex_core->eval(
                  this->mesh->get_coords(vertex),
                  c_data[lid]
                  );
              // Add to matrix
              const auto gid = this->getMap()->getGlobalElement(lid);
              const auto num_lhs = this->sumIntoGlobalValues(
                  gid,
                  Teuchos::tuple<int>(gid),
                  Teuchos::tuple<double>(val.lhs)
                  );
#ifndef NDEBUG
              TEUCHOS_ASSERT_EQUALITY(num_lhs, 1);
#endif
              if (rhs) {
                // add to rhs
                rhs->sumIntoGlobalValue(gid, val.rhs);
              }
            }
          }
        }
      }

      void
      add_domain_boundary_contributions_(
          const std::shared_ptr<Tpetra::Vector<double,int,int>> & rhs
          )
      {
        const auto surfs = this->mesh->boundary_surface_areas();
        for (const auto boundary_core: this->boundary_cores_) {
          for (const auto & subdomain_id: boundary_core->subdomain_ids) {
            const auto verts = this->mesh->get_vertices(subdomain_id);
            for (const auto vert: verts) {
              const auto lid = this->mesh->local_index(vert);
              const auto val = boundary_core->eval(
                  this->mesh->get_coords(vert),
                  surfs[lid]
                  );
              // Add to matrix
              const auto gid = this->getMap()->getGlobalElement(lid);
              const auto num_lhs = this->sumIntoGlobalValues(
                  gid,
                  Teuchos::tuple<int>(gid),
                  Teuchos::tuple<double>(val.lhs)
                  );
#ifndef NDEBUG
              TEUCHOS_ASSERT_EQUALITY(num_lhs, 1);
#endif
              if (rhs) {
                // add to rhs
                rhs->sumIntoGlobalValue(gid, val.rhs);
              }
            }
          }
        }
      }

      // apply dirichlet boundary conditions
      void
      apply_dbcs_(
          const std::shared_ptr<Tpetra::Vector<double,int,int>> & rhs
          )
      {
        for (const auto & bc: this->dbcs_) {
          for (const auto & subdomain_id: bc->subdomain_ids) {
            const auto verts = this->mesh->get_vertices(subdomain_id);
            for (const auto & vertex: verts) {
              // eliminate the row in A
              auto gid = this->mesh->gid(vertex);
              size_t num = this->getNumEntriesInGlobalRow(gid);
              // It shouldn't actually happen that the specified global row
              // does not belong to this graph.
              // TODO find out if why we need this, fix the underlying issue,
              // make this a TEUCHOS_TEST_*
              if (num != Teuchos::OrdinalTraits<size_t>::invalid()) {
                std::vector<int> cols(num);
                std::vector<double> vals(num);
                this->getGlobalRowCopy(gid, cols, vals, num);
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
                this->replaceGlobalValues(gid, cols, vals);
              }
              if (rhs) {
                // set rhs
                const auto coord = this->mesh->get_coords(vertex);
                rhs->replaceGlobalValue(gid, bc->eval(coord));
              }
            }
          }
        }
      }

    public:
      const std::shared_ptr<const nosh::mesh> mesh;

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

#endif // NOSH_FVM_MATRIX_H
