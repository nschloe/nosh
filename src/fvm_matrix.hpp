#ifndef NOSH_FVMMATRIX_H
#define NOSH_FVMMATRIX_H

#include <Teuchos_Tuple.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "matrix.hpp"
#include "matrix_core.hpp"
#include "mesh.hpp"

namespace nosh
{
  class fvm_matrix:
    public matrix
  {
    public:
      fvm_matrix(
          const std::shared_ptr<const nosh::mesh> & _mesh,
          const std::set<std::shared_ptr<const matrix_core>> & matrix_cores,
          const std::set<std::shared_ptr<const nosh::dirichlet_bc>> & _bcs
          ) :
        matrix(_mesh, _bcs)
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        ,fill_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: fvm_matrix::fill_"))
#endif
        // ,matrix_cores_(matrix_cores), TODO
        ,edge_core_associations_(this->associate_edge_cores(matrix_cores))
        ,vertex_core_associations_(this->associate_vertex_cores(matrix_cores))
        {
        }

      virtual
      ~fvm_matrix()
      {};

    protected:
      void
        fill_()
        {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor tm(*fill_time_);
#endif
#ifndef NDEBUG
          TEUCHOS_ASSERT(this->mesh);
#endif
          this->resumeFill();
          this->setAllToScalar(0.0);

          // Add edge contributions
          const std::vector<edge> edges = this->mesh->my_edges();
          const auto edge_data = this->mesh->get_edge_data();
          for (size_t k = 0; k < edges.size(); k++) {
            auto vals = this->edge_core_associations_[k]->edge_contrib(
                this->mesh->get_coords(std::get<0>(edges[k])),
                this->mesh->get_coords(std::get<1>(edges[k])),
                edge_data[k].length,
                edge_data[k].covolume
                );

            const Teuchos::Tuple<int,2> & idx = this->mesh->edge_gids[k];
            for (int i = 0; i < 2; i++) {
              int num = this->sumIntoGlobalValues(
                  idx[i], idx,
                  Teuchos::ArrayView<double>(vals[i])
                  );
#ifndef NDEBUG
              TEUCHOS_ASSERT_EQUALITY(num, 2);
#endif
            }
          }

          // Add vertex contributions
          const auto & control_volumes = this->mesh->control_volumes();
          const auto & owned_nodes = this->mesh->get_owned_nodes();
          auto c_data = control_volumes->getData();
#ifndef NDEBUG
          TEUCHOS_ASSERT_EQUALITY(c_data.size(), owned_nodes.size());
#endif
          for (int k = 0; k < c_data.size(); k++) {
            auto val = this->vertex_core_associations_[k]->vertex_contrib(
                this->mesh->get_coords(owned_nodes[k]),
                c_data[k]
                );
            auto gid = this->getMap()->getGlobalElement(k);
            int num = this->sumIntoGlobalValues(
                gid,
                Teuchos::tuple<int>(gid),
                Teuchos::tuple<double>(val)
                );
#ifndef NDEBUG
            TEUCHOS_ASSERT_EQUALITY(num, 1);
#endif
          }

          this->apply_bcs_();

          this->fillComplete();

          return;
        }
    private:
      // Create a vector core_idx that associated exactly one core with each
      // edge. This core will be the one used for building the matrix at that
      // edge.
      std::vector<std::shared_ptr<const matrix_core>>
        associate_edge_cores(
          const std::set<std::shared_ptr<const matrix_core>> & matrix_cores
          )
      {
        const std::vector<edge> edges = this->mesh->my_edges();
        std::vector<std::shared_ptr<const matrix_core>> core_idx(edges.size());
        for (size_t k = 0; k < edges.size(); k++) {
          const auto x0 = this->mesh->get_coords(std::get<0>(edges[k]));
          const auto x1 = this->mesh->get_coords(std::get<1>(edges[k]));
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
        const auto & owned_nodes = this->mesh->get_owned_nodes();
        std::vector<std::shared_ptr<const matrix_core>> core_idx(owned_nodes.size());
        for (size_t k = 0; k < owned_nodes.size(); k++) {
          const auto x = this->mesh->get_coords(owned_nodes[k]);
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

    private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
      const Teuchos::RCP<Teuchos::Time> fill_time_;
#endif
      std::vector<std::shared_ptr<const matrix_core>> edge_core_associations_;
      std::vector<std::shared_ptr<const matrix_core>> vertex_core_associations_;
  };
} // namespace nosh

#endif // NOSH_FVMMATRIX_H
