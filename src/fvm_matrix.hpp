// @HEADER
//
//    Builder class for the fvm_matrix.
//    Copyright (C) 2015  Nico Schlömer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
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
          const double edge_coefficient,
          const Eigen::Vector3d & edge_midpoint
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
          const auto edge_coefficients = this->mesh->edge_coefficients();
          const vector_fieldType & coords_field =
            this->mesh->get_node_field("coordinates");
          for (size_t k = 0; k < edges.size(); k++) {
            const Eigen::Vector3d & x0 =
              this->mesh->get_node_value(coords_field, std::get<0>(edges[k]));
            const Eigen::Vector3d & x1 =
              this->mesh->get_node_value(coords_field, std::get<1>(edges[k]));
            const Eigen::Vector3d edge_midpoint = 0.5 * (x0 + x1);

            auto vals = edge_contrib(edge_coefficients[k], edge_midpoint);

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
            int num = this->sumIntoLocalValues(
                k,
                Teuchos::tuple<int>(k),
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
