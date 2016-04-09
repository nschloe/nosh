#ifndef NOSH_FVM_HELPERS
#define NOSH_FVM_HELPERS

#include <Teuchos_Tuple.hpp>

#include "helper.hpp"
#include "matrix.hpp"
#include "edge_core.hpp"
#include "vertex_core.hpp"
#include "boundary_core.hpp"
#include "mesh.hpp"

namespace nosh
{
  namespace fvm_helpers
  {
    void
    add_edge_contributions(
        const std::shared_ptr<const nosh::mesh> & mesh,
        const std::set<std::shared_ptr<const edge_core>> & edge_cores,
        const std::shared_ptr<Tpetra::CrsMatrix<double,int,int>> & matrix,
        const std::shared_ptr<Tpetra::Vector<double,int,int>> & rhs
        )
    {
      for (const auto & edge_core: edge_cores) {
        for (const auto & subdomain_id: edge_core->subdomain_ids) {
          const auto edge_data = mesh->get_edge_data();

          // meshset interior edges
          const auto edges = mesh->get_edges(subdomain_id);
          for (const auto edge: edges) {
            const auto verts = mesh->get_vertex_tuple(edge);
            const auto lid = mesh->local_index(edge);
            auto vals = edge_core->eval(
                mesh->get_coords(verts[0]),
                mesh->get_coords(verts[1]),
                edge_data[lid].length,
                edge_data[lid].covolume
                );

            const auto & gids = mesh->edge_gids[lid];
            for (int i = 0; i < 2; i++) {
              if (matrix) {
                // Add to matrix
                const int num_lhs = matrix->sumIntoGlobalValues(
                    gids[i], gids,
                    Teuchos::ArrayView<double>(vals.lhs[i])
                    );
#ifndef NDEBUG
                TEUCHOS_ASSERT_EQUALITY(num_lhs, 2);
#endif
              }
              if (rhs) {
                // Add to rhs
                rhs->sumIntoGlobalValue(gids[i], vals.rhs[i]);
              }
            }
          }

          if (subdomain_id != "everywhere") {
            // meshset boundary edges
            const auto boundary_edges = mesh->get_edges(
                subdomain_id + "_boundary"
                );
            for (const auto edge: boundary_edges) {
              const auto verts = mesh->get_vertex_tuple(edge);

              // check which one of the two verts is in meshset
              int i;
              if (mesh->contains(subdomain_id, {verts[0]})) {
                i = 0;
              } else if (mesh->contains(subdomain_id, {verts[1]})) {
                i = 1;
              } else {
                TEUCHOS_TEST_FOR_EXCEPT_MSG(
                    true,
                    "Neither of the two edge vertices is contained in the subdomain."
                    );
              }

              const auto lid = mesh->local_index(edge);
              auto vals = edge_core->eval(
                  mesh->get_coords(verts[0]),
                  mesh->get_coords(verts[1]),
                  edge_data[lid].length,
                  edge_data[lid].covolume
                  );

              const auto & gids = mesh->edge_gids[lid];
              if (matrix) {
                // Add to matrix
                int num_lhs = matrix->sumIntoGlobalValues(
                    gids[i], gids,
                    Teuchos::ArrayView<double>(vals.lhs[i])
                    );
#ifndef NDEBUG
                TEUCHOS_ASSERT_EQUALITY(num_lhs, 2);
#endif
              }
              if (rhs) {
                // Add to rhs
                rhs->sumIntoGlobalValue(gids[i], vals.rhs[i]);
              }
            }
          }
        }
      }
    }

    void
    add_vertex_contributions(
        const std::shared_ptr<const nosh::mesh> & mesh,
        const std::set<std::shared_ptr<const vertex_core>> & vertex_cores,
        const std::shared_ptr<Tpetra::CrsMatrix<double,int,int>> & matrix,
        const std::shared_ptr<Tpetra::Vector<double,int,int>> & rhs
        )
    {
      const auto & control_volumes = mesh->control_volumes();
      const auto c_data = control_volumes->getData();
      for (const auto & vertex_core: vertex_cores) {
        for (const auto & subdomain_id: vertex_core->subdomain_ids) {
          const auto verts = mesh->get_vertices(subdomain_id);
          for (const auto & vertex: verts) {
            const auto lid = mesh->local_index(vertex);
            const auto val = vertex_core->eval(
                mesh->get_coords(vertex),
                c_data[lid]
                );
            // Add to matrix
            const auto gid = matrix->getMap()->getGlobalElement(lid);
            if (matrix) {
              const auto num_lhs = matrix->sumIntoGlobalValues(
                  gid,
                  Teuchos::tuple<int>(gid),
                  Teuchos::tuple<double>(val.lhs)
                  );
#ifndef NDEBUG
              TEUCHOS_ASSERT_EQUALITY(num_lhs, 1);
#endif
            }
            if (rhs) {
              // add to rhs
              rhs->sumIntoGlobalValue(gid, val.rhs);
            }
          }
        }
      }
    }

    void
    add_domain_boundary_contributions(
        const std::shared_ptr<const nosh::mesh> & mesh,
        const std::set<std::shared_ptr<const boundary_core>> & boundary_cores,
        const std::shared_ptr<Tpetra::CrsMatrix<double,int,int>> & matrix,
        const std::shared_ptr<Tpetra::Vector<double,int,int>> & rhs
        )
    {
      const auto surfs = mesh->boundary_surface_areas();
      for (const auto boundary_core: boundary_cores) {
        for (const auto & subdomain_id: boundary_core->subdomain_ids) {
          const auto verts = mesh->get_vertices(subdomain_id);
          for (const auto vert: verts) {
            const auto lid = mesh->local_index(vert);
            const auto val = boundary_core->eval(
                mesh->get_coords(vert),
                surfs[lid]
                );
            // Add to matrix
            const auto gid = matrix->getMap()->getGlobalElement(lid);
            if (matrix) {
              const auto num_lhs = matrix->sumIntoGlobalValues(
                  gid,
                  Teuchos::tuple<int>(gid),
                  Teuchos::tuple<double>(val.lhs)
                  );
#ifndef NDEBUG
              TEUCHOS_ASSERT_EQUALITY(num_lhs, 1);
#endif
            }
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
    apply_dbcs(
        const std::shared_ptr<const nosh::mesh> & mesh,
        const std::set<std::shared_ptr<const dirichlet_bc>> & dbcs,
        const std::shared_ptr<Tpetra::CrsMatrix<double,int,int>> & matrix,
        const std::shared_ptr<Tpetra::Vector<double,int,int>> & rhs
        )
    {
      for (const auto & bc: dbcs) {
        for (const auto & subdomain_id: bc->subdomain_ids) {
          const auto verts = mesh->get_vertices(subdomain_id);
          for (const auto & vertex: verts) {
            // eliminate the row in A
            auto gid = mesh->gid(vertex);
            if (matrix) {
              size_t num = matrix->getNumEntriesInGlobalRow(gid);
              // It shouldn't actually happen that the specified global row
              // does not belong to this graph.
              // TODO find out if why we need this, fix the underlying issue,
              // make this a TEUCHOS_TEST_*
              if (num != Teuchos::OrdinalTraits<size_t>::invalid()) {
                std::vector<int> cols(num);
                std::vector<double> vals(num);
                matrix->getGlobalRowCopy(gid, cols, vals, num);
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
                matrix->replaceGlobalValues(gid, cols, vals);
              }
            }
            if (rhs) {
              // set rhs
              const auto coord = mesh->get_coords(vertex);
              rhs->replaceGlobalValue(gid, bc->eval(coord));
            }
          }
        }
      }
    }
  } // namespace fvm_helpers
} // namespace nosh

#endif // NOSH_FVM_HELPERS
