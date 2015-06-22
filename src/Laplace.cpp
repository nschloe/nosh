// @HEADER
//
//    Builds the Laplace operator and its variants.
//    Copyright (C) 2012  Nico Schlömer
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
// =============================================================================
// includes
#include "Laplace.hpp"

#include <map>
#include <string>

#include "Mesh.hpp"
#include "ScalarField_Virtual.hpp"

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif

namespace Nosh
{
// =============================================================================
Laplace::
Laplace(
    const std::shared_ptr<const Nosh::Mesh> & _mesh,
    const std::shared_ptr<const Nosh::DirichletBoundaryConditions> & _bcs
    ) :
  LinearOperator(_mesh, _bcs),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  fillTime_(Teuchos::TimeMonitor::getNewTimer("Nosh: Laplace::fill_"))
#endif
{
  this->fill_();
}
// =============================================================================
Laplace::
~Laplace()
{
}
// =============================================================================
void
Laplace::
fill_()
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*fillTime_);
#endif

  this->resumeFill();

  this->setAllToScalar(0.0);

#ifndef NDEBUG
  TEUCHOS_ASSERT(this->mesh);
#endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
  const std::vector<edge> edges = this->mesh->getMyEdges();
  const auto edgeCoefficients = this->mesh->getEdgeCoefficients();

  // Loop over all edges.
  for (size_t k = 0; k < edges.size(); k++) {
    // We'd like to insert the 2x2 matrix
    //
    //     [   alpha, - alpha ]
    //     [ - alpha,   alpha ]
    //
    // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
    // that shares and edge.
    // Do that now, just blockwise for real and imaginary part.
    const double & a = 0.5 * edgeCoefficients[k];
    auto vals = Teuchos::tuple(
      Teuchos::tuple( a, -a),
      Teuchos::tuple(-a,  a)
      );
    const Teuchos::Tuple<int,2> & idx = this->mesh->edgeGids[k];
    for (int i = 0; i < 2; i++) {
      int num = this->sumIntoGlobalValues(idx[i], idx, vals[i]);
#ifndef NDEBUG
      TEUCHOS_ASSERT_EQUALITY(num, 2);
#endif
    }
  }

  this->fillComplete();

  return;
}
// =============================================================================
} // namespace Nosh
