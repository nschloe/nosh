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
#include "ParameterMatrix_Laplace.hpp"

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
namespace ParameterMatrix
{
// =============================================================================
Laplace::
Laplace(
    const std::shared_ptr<const Nosh::Mesh> &mesh,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness
    ) :
  ParameterObject(),
  Tpetra::CrsMatrix<double,int,int>(mesh->buildComplexGraph()),
  mesh_(mesh),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  fillTime_(Teuchos::TimeMonitor::getNewTimer(
               "Nosh: Laplace::fill_")),
#endif
  thickness_(thickness),
  alphaCache_(),
  alphaCacheUpToDate_(false)
{
}
// =============================================================================
Laplace::
~Laplace()
{
}
// =============================================================================
//void
//Laplace::
//apply(const std::map<std::string, double> &params,
//      const Tpetra::Vector<double,int,int> &X,
//      Tpetra::Vector<double,int,int> &Y
//     ) const
//{
//  (void) params;
//  // Rebuild if necessary.
//  if (!matrixCacheUpToDate_) {
//    this->fill_(matrixCache_);
//    matrixCacheUpToDate_ = true;
//  }
//  // This direct application in the cache saves
//  // one matrix copy compared to an explicit
//  // fill(), Apply().
//  TEUCHOS_ASSERT_EQUALITY(0, matrixCache_.Apply(X, Y));
//  return;
//}
// =============================================================================
//void
//Laplace::
//fill(Epetra_FECrsMatrix &matrix,
//     const std::map<std::string, double> &params
//     ) const
//{
//  (void) params;
//  // Cache the construction of the Laplacian.
//  // After all, it's always the same.
//  if (!matrixCacheUpToDate_) {
//    this->fill_(matrixCache_);
//    matrixCacheUpToDate_ = true;
//  }
//
//  matrix = matrixCache_;
//  return;
//}
// =============================================================================
const std::map<std::string, double>
Laplace::
getParameters() const
{
  return std::map<std::string, double>();
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
  TEUCHOS_ASSERT(mesh_);
#endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
#ifndef NDEBUG
  TEUCHOS_ASSERT(thickness_);
#endif

  const std::vector<edge> edges = mesh_->getMyEdges();
  if (!alphaCacheUpToDate_) {
    this->buildAlphaCache_(edges, mesh_->getEdgeCoefficients());
  }

  // Loop over all edges.
  for (std::size_t k = 0; k < edges.size(); k++) {
    // We'd like to insert the 2x2 matrix
    //
    //     [   alpha, - alpha ]
    //     [ - alpha,   alpha ]
    //
    // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
    // that shares and edge.
    // Do that now, just blockwise for real and imaginary part.
    const double & a = alphaCache_[k];
    auto vals = Teuchos::tuple(
      Teuchos::tuple(  a, 0.0,  -a, 0.0),
      Teuchos::tuple(0.0,   a, 0.0,  -a),
      Teuchos::tuple( -a, 0.0,   a, 0.0),
      Teuchos::tuple(0.0,  -a, 0.0,   a)
      );
    const Teuchos::Tuple<int,4> & idx = mesh_->edgeGidsComplex[k];
    for (int i = 0; i < 4; i++) {
      int num = this->sumIntoGlobalValues(idx[i], idx, vals[i]);
#ifndef NDEBUG
      TEUCHOS_ASSERT_EQUALITY(num, 4);
#endif
    }
  }

  this->fillComplete();

  return;
}
// =============================================================================
void
Laplace::
buildAlphaCache_(
    const std::vector<edge> & edges,
    const std::vector<double> & edgeCoefficients
    ) const
{
  // This routine serves the one and only purpose of caching the thickness
  // average. The cache is used in every call to this->fill().  This is
  // somewhat problematic since the map of V is principally not known here.
  // Also, it is typically a nonoverlapping map whereas some edges do sit on a
  // processor boundary, so actually the values of V are needed in an
  // overlapping map.
  // Fair enough. Let's distribute the vales of V to an overlapping map here.
  alphaCache_ = std::vector<double>(edges.size());

  std::map<std::string, double> dummy;
  const auto thicknessValues = thickness_->getV(dummy);

  auto overlapMap = mesh_->getOverlapMap();
  // We need to make sure that thicknessValues are distributed on the overlap
  // map.
  // Make sure to use Import here instead of Export as the vector that we want
  // to build is overlapping, "larger". If the "smaller", non-overlapping
  // vector is exported, only the values on the overlap would only be set on
  // one processor.
  Tpetra::Vector<double,int,int> thicknessOverlap(Teuchos::rcp(overlapMap));
  Tpetra::Import<int,int> importer(
      thicknessValues.getMap(),
      Teuchos::rcp(overlapMap)
      );
  thicknessOverlap.doImport(thicknessValues, importer, Tpetra::INSERT);

  Teuchos::ArrayRCP<const double> tData = thicknessOverlap.getData();

  int gid0, gid1;
  int lid0, lid1;
  for (std::size_t k = 0; k < edges.size(); k++) {
    // Get the ID of the edge endpoints in the map of getV(). Well...
    gid0 = mesh_->gid(std::get<0>(edges[k]));
    lid0 = overlapMap->getLocalElement(gid0);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid0 < 0,
        "The global index " << gid0
        << " does not seem to be present on this node."
       );
#endif
    gid1 = mesh_->gid(std::get<1>(edges[k]));
    lid1 = overlapMap->getLocalElement(gid1);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid1 < 0,
        "The global index " << gid1
        << " does not seem to be present on this node."
       );
#endif
    // Update cache.
    alphaCache_[k] = edgeCoefficients[k]
      * 0.5 * (tData[lid0] + tData[lid1]);
  }

  alphaCacheUpToDate_ = true;

  return;
}
// =============================================================================
} // namespace ParameterMatrix
} // namespace Nosh
