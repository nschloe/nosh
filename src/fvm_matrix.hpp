#ifndef NOSH_FVMMATRIX_H
#define NOSH_FVMMATRIX_H

#include <Teuchos_Tuple.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "matrix.hpp"
#include "mesh.hpp"

namespace nosh
{
  class fvm_matrix:
    public matrix
  {
    public:
      fvm_matrix(
          const std::shared_ptr<const nosh::mesh> & _mesh,
          const std::set<std::shared_ptr<const nosh::dirichlet_bc>> & _bcs
          ) :
        matrix(_mesh, _bcs)
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        ,fill_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: fvm_matrix::fill_"))
#endif
        {
        }

      virtual
      ~fvm_matrix()
      {};

    protected:

      virtual
      std::vector<std::vector<double>>
      edge_contrib(
          const Eigen::Vector3d & x0,
          const Eigen::Vector3d & x1,
          const double edge_length,
          const double covolume
          ) const = 0;

      virtual
      double
      vertex_contrib(
          const double control_volume
          ) const = 0;

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
            auto vals = edge_contrib(
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
          auto c_data = control_volumes->getData();
          for (int k = 0; k < c_data.size(); k++) {
            auto val = vertex_contrib(c_data[k]);
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
#ifdef NOSH_TEUCHOS_TIME_MONITOR
      const Teuchos::RCP<Teuchos::Time> fill_time_;
#endif
  };
} // namespace nosh

#endif // NOSH_FVMMATRIX_H
