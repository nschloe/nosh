// @HEADER
//
//    Builder class for the EdgeMatrix.
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
#ifndef NOSH_EDGEOPERATOR_H
#define NOSH_EDGEOPERATOR_H

#include <Teuchos_Tuple.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "LinearOperator.hpp"
#include "Mesh.hpp"

namespace Nosh
{
  class EdgeMatrix:
    public LinearOperator
  {
    public:
      EdgeMatrix(
          const std::shared_ptr<const Nosh::Mesh> & _mesh,
          const std::set<std::shared_ptr<const Nosh::DirichletBC>> & _bcs
          ) :
        LinearOperator(_mesh, _bcs)
#ifdef NOSH_TEUCHOS_TIME_MONITOR
        ,fillTime_(Teuchos::TimeMonitor::getNewTimer("Nosh: EdgeMatrix::fill_"))
#endif
        {
        }

      virtual
      ~EdgeMatrix()
      {};

    protected:

      virtual
      std::vector<std::vector<double>>
      edgeContrib(
          const double edgeCoefficient,
          const double controlVolume0,
          const double controlVolume1
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
          const std::vector<edge> edges = this->mesh->getMyEdges();
          const auto edgeCoefficients = this->mesh->getEdgeCoefficients();
          for (size_t k = 0; k < edges.size(); k++) {
            const double & a = edgeCoefficients[k];
            auto vals = edgeContrib(
                a, 0.0, 0.0 // TODO provide correct values
                );
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

#endif // NOSH_EDGEOPERATOR_H
