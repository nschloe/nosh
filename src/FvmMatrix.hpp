// @HEADER
//
//    Builder class for the FvmMatrix.
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

#include "LinearOperator.hpp"
#include "Mesh.hpp"

namespace Nosh
{
  class FvmMatrix:
    public LinearOperator
  {
    public:
      FvmMatrix(
          const std::shared_ptr<const Nosh::Mesh> & _mesh,
          const std::set<std::shared_ptr<const Nosh::DirichletBC>> & _bcs
          ) :
        LinearOperator(_mesh, _bcs)
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        ,fillTime_(Teuchos::TimeMonitor::getNewTimer("Nosh: FvmMatrix::fill_"))
#endif
        {
        }

      virtual
      ~FvmMatrix()
      {};

    protected:

      virtual
      std::vector<std::vector<double>>
      edgeContrib(
          const double edgeCoefficient,
          const Eigen::Vector3d & edgeMidpoint
          ) const = 0;

      virtual
      double
      vertexContrib(
          const double controlVolume
          ) const = 0;

      void
        fill_()
        {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor tm(*fillTime_);
#endif
#ifndef NDEBUG
          TEUCHOS_ASSERT(this->mesh);
#endif
          this->resumeFill();
          this->setAllToScalar(0.0);

          // Add edge contributions
          const std::vector<edge> edges = this->mesh->getMyEdges();
          const auto edgeCoefficients = this->mesh->getEdgeCoefficients();
          const VectorFieldType & coordsField =
            this->mesh->getNodeField("coordinates");
          for (size_t k = 0; k < edges.size(); k++) {
            const Eigen::Vector3d & x0 =
              this->mesh->getNodeValue(coordsField, std::get<0>(edges[k]));
            const Eigen::Vector3d & x1 =
              this->mesh->getNodeValue(coordsField, std::get<1>(edges[k]));
            const Eigen::Vector3d edgeMidpoint = 0.5 * (x0 + x1);

            auto vals = edgeContrib(edgeCoefficients[k], edgeMidpoint);

            const Teuchos::Tuple<int,2> & idx = this->mesh->edgeGids[k];
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
          const auto & controlVolumes = this->mesh->getControlVolumes();
          auto cData = controlVolumes->getData();
          for (int k = 0; k < cData.size(); k++) {
            auto val = vertexContrib(cData[k]);
            int num = this->sumIntoLocalValues(
                k,
                Teuchos::tuple<int>(k),
                Teuchos::tuple<double>(val)
                );
#ifndef NDEBUG
            TEUCHOS_ASSERT_EQUALITY(num, 1);
#endif
          }

          this->applyBcs_();

          this->fillComplete();

          return;
        }

    private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
      const Teuchos::RCP<Teuchos::Time> fillTime_;
#endif
  };
} // namespace Nosh

#endif // NOSH_FVMMATRIX_H
