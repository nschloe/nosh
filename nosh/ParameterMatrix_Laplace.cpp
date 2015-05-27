// @HEADER
//
//    Builds the Laplace operator and its variants.
//    Copyright (C) 2012  Nico Schl\"omer
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
#include "nosh/ParameterMatrix_Laplace.hpp"

#include <map>
#include <string>

#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Virtual.hpp"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>

#include <ml_epetra_preconditioner.h>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif

namespace Nosh
{
namespace ParameterMatrix
{
// =============================================================================
Laplace::
Laplace(const Teuchos::RCP<const Nosh::StkMesh> &mesh,
        const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness
        ) :
  Virtual(mesh),
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
//      const Epetra_Vector &X,
//      Epetra_Vector &Y
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
//// =============================================================================
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
  // Zero-out the matrix.
  TEUCHOS_ASSERT_EQUALITY(0, this->PutScalar(0.0));

#ifndef NDEBUG
  TEUCHOS_ASSERT(!mesh_.is_null());
#endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
#ifndef NDEBUG
  TEUCHOS_ASSERT(!thickness_.is_null());
#endif

  const std::vector<edge> edges = mesh_->getEdgeNodes();
  if (!alphaCacheUpToDate_)
    this->buildAlphaCache_(edges, mesh_->getEdgeCoefficients());

  // Loop over all edges.
  for (auto k = 0; k < edges.size(); k++) {
    // We'd like to insert the 2x2 matrix
    //
    //     [   alpha, - alpha ]
    //     [ - alpha,   alpha ]
    //
    // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
    // that shares and edge.
    // Do that now, just blockwise for real and imaginary part.
    const double & a = alphaCache_[k];
    double ain [] = {
        a, 0.0,  -a, 0.0,
      0.0,   a, 0.0,  -a,
       -a, 0.0,   a, 0.0,
      0.0,  -a, 0.0,   a
    };
    int ierr = this->SumIntoGlobalValues(
        mesh_->globalIndexCache[k],
        Epetra_SerialDenseMatrix(View, ain, 4, 4 ,4)
        );
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(0, ierr);
#endif
  }

  // calls FillComplete by default
  TEUCHOS_ASSERT_EQUALITY(0, this->GlobalAssemble());
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
  alphaCache_ = std::vector<double>(edges.size());

  std::map<std::string, double> dummy;
  const Epetra_Vector thicknessValues = thickness_->getV(dummy);

  int gid0, gid1;
  int lid0, lid1;
  for (auto k = 0; k < edges.size(); k++) {
    gid0 = mesh_->gid(std::get<0>(edges[k]));
    lid0 = mesh_->getNodesOverlapMap()->LID(gid0);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid0 < 0,
        "The global index " << gid0
        << " does not seem to be present on this node."
        );
#endif
    gid1 = mesh_->gid(std::get<1>(edges[k]));
    lid1 = mesh_->getNodesOverlapMap()->LID(gid1);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid1 < 0,
        "The global index " << gid1
        << " does not seem to be present on this node."
        );
#endif
    alphaCache_[k] = edgeCoefficients[k]
      * 0.5 * (thicknessValues[lid0] + thicknessValues[lid1]);
  }

  alphaCacheUpToDate_ = true;

  return;
}
// =============================================================================
} // namespace ParameterMatrix
} // namespace Nosh
