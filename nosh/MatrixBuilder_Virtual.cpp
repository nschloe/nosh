// @HEADER
//
//    Virtual class for matrix constructors.
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

#include "nosh/MatrixBuilder_Virtual.hpp"

#include "nosh/StkMesh.hpp"

namespace Nosh
{
namespace MatrixBuilder
{
// ============================================================================
Virtual::
Virtual(const Teuchos::RCP<const Nosh::StkMesh> &mesh):
  mesh_(mesh),
  globalIndexCache_(this->buildGlobalIndexCache_()),
  graph_(this->buildGraph_())
{
}
// ============================================================================
Virtual::
~Virtual()
{
}
// =============================================================================
const Epetra_Comm &
Virtual::
getComm() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(!mesh_.is_null());
#endif
  return mesh_->getComm();
}
// =============================================================================
const Epetra_FECrsGraph &
Virtual::
getGraph() const
{
  return graph_;
}
// =============================================================================
const Epetra_FECrsGraph
Virtual::
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

  const Teuchos::Array<edge> edges = mesh_->getEdgeNodes();

  // Loop over all edges and put entries wherever two nodes are connected.
  for (auto k = 0; k < edges.size(); k++) {
    int ierr = graph.InsertGlobalIndices(
        4, globalIndexCache_[k].Values(),
        4, globalIndexCache_[k].Values()
        );
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(0, ierr);
#endif
  }

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  int ierr = graph.GlobalAssemble(noMap, noMap);
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(0, ierr);
#endif

  return graph;
}
// =============================================================================
const Teuchos::ArrayRCP<Epetra_IntSerialDenseVector>
Virtual::
buildGlobalIndexCache_() const
{
  const Teuchos::Array<edge> edges = mesh_->getEdgeNodes();

  Teuchos::ArrayRCP<Epetra_IntSerialDenseVector> globalIndexCache(edges.size());

  Teuchos::Tuple<int, 2> gid;
  for (auto k = 0; k < edges.size(); k++) {
    gid[0] = mesh_->gid(std::get<0>(edges[k]));
    gid[1] = mesh_->gid(std::get<1>(edges[k]));

    globalIndexCache_[k] = Epetra_IntSerialDenseVector(4);
    globalIndexCache_[k][0] = 2 * gid[0];
    globalIndexCache_[k][1] = 2 * gid[0] + 1;
    globalIndexCache_[k][2] = 2 * gid[1];
    globalIndexCache_[k][3] = 2 * gid[1] + 1;
  }

  return globalIndexCache;
}
// ============================================================================
} // namespace MatrixBuilder
} // namespace Nosh
