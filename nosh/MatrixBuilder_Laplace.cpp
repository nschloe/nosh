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
#include "nosh/MatrixBuilder_Laplace.hpp"

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
namespace MatrixBuilder
{
// =============================================================================
Laplace::
Laplace(const Teuchos::RCP<const Nosh::StkMesh> &mesh,
        const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness
     ) :
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  fillTime_(Teuchos::TimeMonitor::getNewTimer(
               "Nosh: Laplace::fill_")),
#endif
  mesh_(mesh),
  thickness_(thickness),
  globalIndexCache_(),
  globalIndexCacheUpToDate_(false),
  graph_(this->buildGraph_()), // build the graph immediately
  matrixCache_(Copy, graph_),
  matrixCacheUpToDate_(false),
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
const Epetra_Comm &
Laplace::
getComm() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(!mesh_.is_null());
#endif
  return mesh_->getComm();
}
// =============================================================================
const Epetra_FECrsGraph &
Laplace::
getGraph() const
{
  return graph_;
}
// =============================================================================
void
Laplace::
apply(const std::map<std::string, double> &params,
      const Epetra_Vector &X,
      Epetra_Vector &Y
     ) const
{
  (void) params;
  // Rebuild if necessary.
  if (!matrixCacheUpToDate_) {
    this->fill_(matrixCache_);
    matrixCacheUpToDate_ = true;
  }
  // This direct application in the cache saves
  // one matrix copy compared to an explicit
  // fill(), Apply().
  TEUCHOS_ASSERT_EQUALITY(0, matrixCache_.Apply(X, Y));
  return;
}
// =============================================================================
void
Laplace::
applyDKDp(const std::map<std::string, double> &params,
          const std::string & paramName,
          const Epetra_Vector &X,
          Epetra_Vector &Y
       ) const
{
  (void) params;
  (void) paramName;
  (void) X;
  (void) Y;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                              "No parameters for the Laplace operator.");
  return;
}
// =============================================================================
void
Laplace::
fill(Epetra_FECrsMatrix &matrix,
     const std::map<std::string, double> &params
  ) const
{
  (void) params;
  // Cache the construction of the Laplacian.
  // After all, it's always the same.
  if (!matrixCacheUpToDate_) {
    this->fill_(matrixCache_);
    matrixCacheUpToDate_ = true;
  }

  matrix = matrixCache_;
  return;
}
// =============================================================================
const std::map<std::string, double>
Laplace::
getInitialParameters() const
{
  return std::map<std::string, double>();
}
// =============================================================================
const Epetra_FECrsGraph
Laplace::
buildGraph_() const
{
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from
  // the ownedNodes map, and the overlapping one from the
  // overlapNodes.
  // Let's illustrate the implications with the example of the matrix
  //   [ 2 1   ]
  //   [ 1 2 1 ]
  //   [   1 2 ].
  // Suppose subdomain 1 consists of node 1, subdomain 2 of node 3,
  // and node 2 forms the boundary between them.
  // For two processes, if process 1 owns nodes 1 and 2, the matrix
  // will be split as
  //   [ 2 1   ]   [       ]
  //   [ 1 2 1 ] + [       ]
  //   [       ]   [   1 2 ].
  // The vectors always need to have a unique map (otherwise, norms
  // cannot be computed by Epetra), so let's assume they have the
  // map ([1, 2], [3]).
  // The communucation for a matrix-vector multiplication Ax=y
  // needs to be:
  //
  //   1. Communicate x(3) to process 1.
  //   2. Communicate x(2) to process 2.
  //   3. Compute.
  //
  // If the matrix is split up like
  //   [ 2 1   ]   [       ]
  //   [ 1 1   ] + [   1 1 ]
  //   [       ]   [   1 2 ]
  // (like the overlap map suggests), then any Ax=y comes down to:
  //
  //   1. Communicate x(2) to process 2.
  //   2. Compute.
  //   3. Communicate (part of) y(2) to process 1.
  //
  // In the general case, assuming that the number of nodes adjacent
  // to a boundary (on one side) are approximately the number of
  // nodes on that boundary, there is not much difference in
  // communication between the patterns.
  // What does differ, though, is the workload on the processes
  // during the computation phase: Process 1 that owns the whole
  // boundary, has to compute more than process 2.
  // Notice, however, that the total number of computations is
  // lower in scenario 1 (7 vs. 8 FLOPs); the same is true for
  // storage.
  // Hence, it comes down to the question whether or not the
  // mesh generator provided a fair share of the boundary nodes.
  // If yes, then scenario 1 will yield approximately even
  // computation times; if not, then scenario 2 will guarantee
  // equal computation times at the cost of higher total
  // storage and computation needs.
  //
  // Remark:
  // This matrix will later be fed into ML. ML has certain restrictions as to
  // what maps can be used. One of those is that RowMatrixRowMap() and
  // OperatorRangeMap must be the same, and, if the matrix is square,
  // OperatorRangeMap and OperatorDomainMap must coincide too.
  //
#ifndef NDEBUG
  TEUCHOS_ASSERT(!mesh_.is_null());
#endif
  const Epetra_Map &noMap = *mesh_->getComplexNonOverlapMap();
  Epetra_FECrsGraph graph(Copy, noMap, 0);

  const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity, 2> > edges =
    mesh_->getEdgeNodes();
  if (!globalIndexCacheUpToDate_)
    this->buildGlobalIndexCache_(edges);

  // Loop over all edges and put entries wherever two nodes are connected.
  for (Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity, 2> >::size_type k = 0;
       k < edges.size();
       k++)
    TEUCHOS_ASSERT_EQUALITY(
        0,
        graph.InsertGlobalIndices(
          4, globalIndexCache_[k].Values(), 4, globalIndexCache_[k].Values()
         )
       );

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  TEUCHOS_ASSERT_EQUALITY(0, graph.GlobalAssemble(noMap, noMap));

  return graph;
}
// =============================================================================
void
Laplace::
fill_(Epetra_FECrsMatrix &matrix) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*fillTime_);
#endif
  // Zero-out the matrix.
  TEUCHOS_ASSERT_EQUALITY(0, matrix.PutScalar(0.0));

#ifndef NDEBUG
  TEUCHOS_ASSERT(!mesh_.is_null());
#endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
#ifndef NDEBUG
  TEUCHOS_ASSERT(!thickness_.is_null());
#endif

  const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity, 2> > edges =
    mesh_->getEdgeNodes();
  if (!globalIndexCacheUpToDate_)
    this->buildGlobalIndexCache_(edges);
  if (!alphaCacheUpToDate_)
    this->buildAlphaCache_(edges, mesh_->getEdgeCoefficients());

  Epetra_SerialDenseMatrix A(4, 4);
  // Loop over all edges.
  for (Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity, 2> >::size_type k = 0;
       k < edges.size();
       k++) {
    // We'd like to insert the 2x2 matrix
    //
    //     [   alpha, - alpha ]
    //     [ - alpha,   alpha ]
    //
    // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
    // that shares and edge.
    // Do that now, just blockwise for real and imaginary part.
    const double & a = alphaCache_[k];
    A(0, 0) =  a;
    A(0, 1) =  0.0;
    A(0, 2) = -a;
    A(0, 3) =  0.0;
    A(1, 0) =  0.0;
    A(1, 1) =  a;
    A(1, 2) =  0.0;
    A(1, 3) = -a;
    A(2, 0) = -a;
    A(2, 1) =  0.0;
    A(2, 2) =  a;
    A(2, 3) =  0.0;
    A(3, 0) =  0.0;
    A(3, 1) = -a;
    A(3, 2) =  0.0;
    A(3, 3) =  a;
    TEUCHOS_ASSERT_EQUALITY(
        0,
        matrix.SumIntoGlobalValues(globalIndexCache_[k], A)
        );
    // -------------------------------------------------------------------
  }

  // calls FillComplete by default
  TEUCHOS_ASSERT_EQUALITY(0, matrix.GlobalAssemble());
  return;
}
// =============================================================================
void
Laplace::
buildGlobalIndexCache_(const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity, 2> > &edges) const
{
  globalIndexCache_ =
    Teuchos::ArrayRCP<Epetra_IntSerialDenseVector>(edges.size());

  Teuchos::Tuple<int, 2> gid;
  for (Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity, 2> >::size_type k = 0;
       k < edges.size();
       k++) {
    gid[0] = edges[k][0]->identifier() - 1;
    gid[1] = edges[k][1]->identifier() - 1;

    globalIndexCache_[k] = Epetra_IntSerialDenseVector(4);
    globalIndexCache_[k][0] = 2*gid[0];
    globalIndexCache_[k][1] = 2*gid[0]+1;
    globalIndexCache_[k][2] = 2*gid[1];
    globalIndexCache_[k][3] = 2*gid[1]+1;
  }

  globalIndexCacheUpToDate_ = true;

  return;
}
// =============================================================================
void
Laplace::
buildAlphaCache_(const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity, 2> > & edges,
                  const Teuchos::ArrayRCP<const double> &edgeCoefficients
               ) const
{
  alphaCache_ = Teuchos::ArrayRCP<double>(edges.size());

  std::map<std::string, double> dummy;
  const Epetra_Vector thicknessValues = thickness_->getV(dummy);

  Teuchos::Tuple<int, 2> gid;
  Teuchos::Tuple<int, 2> lid;
  for (Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity, 2> >::size_type k = 0;
       k < edges.size();
       k++) {
    gid[0] = edges[k][0]->identifier() - 1;
    lid[0] = mesh_->getNodesOverlapMap()->LID(gid[0]);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(lid[0] < 0,
                                 "The global index " << gid[0]
                                 << " does not seem to be present on this node.");
#endif
    gid[1] = edges[k][1]->identifier() - 1;
    lid[1] = mesh_->getNodesOverlapMap()->LID(gid[1]);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(lid[1] < 0,
                                "The global index " << gid[1]
                                << " does not seem to be present on this node.");
#endif
    alphaCache_[k] = edgeCoefficients[k]
                     * 0.5 * (thicknessValues[lid[0]] + thicknessValues[lid[1]]);
  }

  alphaCacheUpToDate_ = true;

  return;
}
// =============================================================================
} // namespace MatrixBuilder
} // namespace Nosh
