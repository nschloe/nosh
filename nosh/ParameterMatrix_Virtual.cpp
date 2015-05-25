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

#include "nosh/ParameterMatrix_Virtual.hpp"

#include "nosh/StkMesh.hpp"

namespace Nosh
{
namespace ParameterMatrix
{
// ============================================================================
Virtual::
Virtual(const Teuchos::RCP<const Nosh::StkMesh> &mesh):
  Epetra_FECrsMatrix(Copy, this->buildGraph_(*mesh)),
  mesh_(mesh),
  buildParameters_()
{
}
// ============================================================================
Virtual::
~Virtual()
{
}
// ============================================================================
void
Virtual::
refill(const std::map<std::string, double> &params)
{
  // Cache the construction of the matrix.
  // This is useful because in the continuation context, the matrix is called a
  // number of times with the same arguments (in computeF, getJacobian(), and
  // getPreconditioner().
  bool needsRefill;
  if (buildParameters_.empty()) {
    needsRefill = true;
  } else {
    needsRefill = false;
    for (auto it = buildParameters_.begin();
         it != buildParameters_.end();
         ++it) {
      // Check if it->first is in params at all and if their values are equal.
      std::map<std::string, double>::const_iterator it2 =
        params.find(it->first);
      TEUCHOS_ASSERT(it2 != params.end());
      if (it2->second != it->second) {
        needsRefill = true;
        break;
      }
    }
  }

  if (needsRefill) {
    this->refill_(params);
    // Reset build parameters.
    for (auto it = buildParameters_.begin();
         it != buildParameters_.end();
         ++it) {
      std::map<std::string, double>::const_iterator it2 =
        params.find(it->first);
      TEUCHOS_ASSERT(it2 != params.end());
      it->second = it2->second;
    }
  }

  return;
}
// =============================================================================
const Epetra_FECrsGraph
Virtual::
buildGraph_(const Nosh::StkMesh & mesh)
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
  const Teuchos::RCP<const Epetra_Map> noMap = mesh.getComplexNonOverlapMap();
#ifndef NDEBUG
    TEUCHOS_ASSERT(!noMap.is_null());
#endif
  Epetra_FECrsGraph graph(Copy, *noMap, 0);

  const Teuchos::Array<edge> edges = mesh.getEdgeNodes();

  // Loop over all edges and put entries wherever two nodes are connected.
  for (auto k = 0; k < edges.size(); k++) {
    int ierr = graph.InsertGlobalIndices(
        4, mesh.globalIndexCache[k].Values(),
        4, mesh.globalIndexCache[k].Values()
        );
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(0, ierr);
#endif
  }

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  int ierr = graph.GlobalAssemble(*noMap, *noMap);
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(0, ierr);
#endif

  return graph;
}
// ============================================================================
} // namespace ParameterMatrix
} // namespace Nosh
