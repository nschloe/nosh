// @HEADER
//
//    Mesh class with compatibility to stk_mesh.
//    Copyright (C) 2010--2012  Nico Schlömer
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
#include "MeshTri.hpp"

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <stk_mesh/base/MetaData.hpp>

namespace Nosh
{
// =============================================================================
MeshTri::
MeshTri(
    const std::shared_ptr<const Teuchos::Comm<int>> & _comm,
    const std::shared_ptr<stk::io::StkMeshIoBroker> & broker
    ) :
  Mesh(_comm, broker),
  controlVolumes_(this->computeControlVolumes_()),
  edgeCoefficients_(this->computeEdgeCoefficients_()),
  boundaryNodes_(this->computeBoundaryNodes_())
{
}
// =============================================================================
MeshTri::
~MeshTri()
{
}
// =============================================================================
std::vector<double>
MeshTri::
computeEdgeCoefficients_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*computeEdgeCoefficientsTime_);
#endif

  std::vector<stk::mesh::Entity> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  auto numEdges = edgeData_.edgeNodes.size();

  std::vector<double> edgeCoefficients(numEdges);

  const VectorFieldType & coordsField = getNodeField("coordinates");

  // Calculate the contributions edge by edge.
  for (unsigned int k = 0; k < numCells; k++) {
    // Get edge coordinates.
    size_t numLocalEdges = edgeData_.cellEdges[k].size();
    std::vector<Eigen::Vector3d> localEdgeCoords(numLocalEdges);
    for (size_t i = 0; i < numLocalEdges; i++) {
      const edge & e = edgeData_.edgeNodes[edgeData_.cellEdges[k][i]];
      localEdgeCoords[i] =
        this->getNodeValue(coordsField, std::get<1>(e))
        - this->getNodeValue(coordsField, std::get<0>(e));
    }

    Eigen::VectorXd edgeCoeffs =
      getEdgeCoefficientsNumerically_(localEdgeCoords);

    // Fill the edge coefficients into the vector.
    for (size_t i = 0; i < numLocalEdges; i++) {
      edgeCoefficients[edgeData_.cellEdges[k][i]] += edgeCoeffs[i];
    }
  }

  return edgeCoefficients;
}
// =============================================================================
Eigen::VectorXd
MeshTri::
getEdgeCoefficientsNumerically_(
  const std::vector<Eigen::Vector3d> edges
  ) const
{
  size_t numEdges = edges.size();

  // Build an equation system for the edge coefficients alpha_k.
  // They fulfill
  //
  //    |simplex| * <u,v> = \sum_{edges e_i} alpha_i <u,e_i> <e_i,v>
  //
  // for any pair of vectors u, v in the plane of the triangle.
  //
  const double vol = 0.5 * (edges[0].cross(edges[1])).norm();

  Eigen::MatrixXd A(numEdges, numEdges);
  Eigen::VectorXd rhs(numEdges);

  // Build the equation system:
  // The equation
  //
  //    |simplex| ||u||^2 = \sum_i \alpha_i <u,e_i> <e_i,u>
  //
  // has to hold for all vectors u in the plane spanned by the edges,
  // particularly by the edges themselves.
  //
  // Only fill the upper part of the Hermitian matrix.
  //
  for (size_t i = 0; i < numEdges; i++) {
    double alpha = edges[i].dot(edges[i]);
    rhs(i) = vol * alpha;
    A(i,i) = alpha * alpha;
    for (size_t j = i+1; j < numEdges; j++) {
      A(i, j) = edges[i].dot(edges[j]) * edges[j].dot(edges[i]);
      A(j, i) = A(i, j);
    }
  }

  // Solve the equation system for the alpha_i.  The system is symmetric and,
  // if the simplex is not degenerate, positive definite.
  //return A.ldlt().solve(rhs);
  return A.fullPivLu().solve(rhs);
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
MeshTri::
computeControlVolumes_() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(nodesMap_);
  TEUCHOS_ASSERT(nodesOverlapMap_);
  TEUCHOS_ASSERT(ioBroker_);
#endif

  // Create temporaries to hold the overlap values for control volumes.
  Tpetra::Vector<double,int,int> cvOverlap(Teuchos::rcp(nodesOverlapMap_));

  this->computeControlVolumesT_(cvOverlap);

  // Export control volumes to a non-overlapping map, and sum the entries.
  Teuchos::RCP<const Tpetra::Export<int,int>> exporter = Tpetra::createExport(
      Teuchos::rcp(nodesOverlapMap_),
      Teuchos::rcp(nodesMap_)
      );

  auto controlVolumes = std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(nodesMap_)
      );
  controlVolumes->doExport(cvOverlap, *exporter, Tpetra::ADD);

  return controlVolumes;
}
// =============================================================================
void
MeshTri::
computeControlVolumesT_(Tpetra::Vector<double,int,int> & cvOverlap) const
{
  std::vector<stk::mesh::Entity> cells = this->getOwnedCells();
  size_t numCells = cells.size();

  Teuchos::ArrayRCP<double> cvData = cvOverlap.getDataNonConst();

  const VectorFieldType & coordsField = getNodeField("coordinates");

  // Calculate the contributions to the finite volumes cell by cell.
  for (size_t k = 0; k < numCells; k++) {
    const stk::mesh::Entity * localNodes =
      ioBroker_->bulk_data().begin_nodes(cells[k]);
    unsigned int numLocalNodes = ioBroker_->bulk_data().num_nodes(cells[k]);

#ifndef NDEBUG
    // Confirm that we always have the same simplices.
    TEUCHOS_ASSERT_EQUALITY(numLocalNodes, 3);
#endif

    // Fetch the nodal positions into 'localNodes'.
    //const std::vector<Eigen::Vector3d> localNodeCoords =
    //this->getNodeCoordinates_(localNodes);
    std::vector<Eigen::Vector3d> localNodeCoords(numLocalNodes);
    for (unsigned int i = 0; i < numLocalNodes; i++) {
      localNodeCoords[i] = this->getNodeValue(coordsField, localNodes[i]);
    }

    // compute the circumcenter of the cell
    const Eigen::Vector3d cc = computeTriangleCircumcenter_(localNodeCoords);

    // Iterate over the edges.
    // As true edge entities are not available here, loop over all pairs of
    // local nodes.
    for (unsigned int e0 = 0; e0 < numLocalNodes; e0++) {
      const Eigen::Vector3d &x0 = localNodeCoords[e0];
      const int gid0 = ioBroker_->bulk_data().identifier(localNodes[e0]) - 1;
      const int lid0 = nodesOverlapMap_->getLocalElement(gid0);
#ifndef NDEBUG
      TEUCHOS_ASSERT_INEQUALITY(lid0, >=, 0);
#endif
      for (unsigned int e1 = e0+1; e1 < numLocalNodes; e1++) {
        const Eigen::Vector3d &x1 = localNodeCoords[e1];
        const int gid1 = ioBroker_->bulk_data().identifier(localNodes[e1]) - 1;
        const int lid1 = nodesOverlapMap_->getLocalElement(gid1);
#ifndef NDEBUG
        TEUCHOS_ASSERT_INEQUALITY(lid1, >=, 0);
#endif
        // Get the other node.
        const unsigned int other = this->getOtherIndex_(e0, e1);

        double edgeLength = (x1-x0).norm();

        // Compute the (n-1)-dimensional covolume.
        double covolume;
        const Eigen::Vector3d &other0 = localNodeCoords[other];
        covolume = this->computeCovolume2d_(cc, x0, x1, other0);
        // The problem with counting the average thickness in 2D is the
        // following.  Ideally, one would want to loop over all edges, add the
        // midpoint value of the thickness to both of the edge end points, and
        // eventually loop over all endpoints and divide by the number of edges
        // (connections) they have with neighboring nodes).
        // Unfortunately, this is impossible now b/c there's no edge generation
        // for shells in Trilinos yet (2011-04-15).
        // As a workaround, one could loop over all cells, and then all pairs
        // of nodes to retrieve the edges. In 2D, almost all of the edges would
        // be counted twice this way as they belong to two cells. This is true
        // for all but the boundary edges. Again, it is difficult (impossible?)
        // to know what the boundary edges are, and hence which values to
        // divide by 2. Dividing them all by two would result in an
        // artificially lower thickness near the boundaries. This is not what
        // we want.

        // Compute the contributions to the finite volumes of the adjacent
        // edges.
        double pyramidVolume = 0.5*edgeLength * covolume / 2;
        cvData[lid0] += pyramidVolume;
        cvData[lid1] += pyramidVolume;
      }
    }
  }

  return;
}
// =============================================================================
double
MeshTri::
computeCovolume2d_(
    const Eigen::Vector3d &cc,
    const Eigen::Vector3d &x0,
    const Eigen::Vector3d &x1,
    const Eigen::Vector3d &other0
    ) const
{
  // edge midpoint
  Eigen::Vector3d mp = 0.5 * (x0 + x1);

  double coedgeLength = (mp - cc).norm();

  // The only difficulty here is to determine whether the length of coedge is
  // to be taken positive or negative.
  // To this end, make sure that the order (x0, cc, mp) is of the same
  // orientation as (x0, other0, mp).
  Eigen::Vector3d cellNormal = (other0 - x0).cross(mp - x0);
  Eigen::Vector3d ccNormal = (cc - x0).cross(mp - x0);

  // copysign takes the absolute value of the first argument and the sign of
  // the second.
  return copysign(coedgeLength, ccNormal.dot(cellNormal));
}
// =============================================================================
unsigned int
MeshTri::
getOtherIndex_(unsigned int e0, unsigned int e1) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_INEQUALITY(e0, !=, e1);
#endif

  // Get the index in [0,1,2] which is not e0, e1.
  if (0 != e0 && 0 != e1)
    return 0;
  else if (1 != e0 && 1 != e1)
    return 1;
  else if (2 != e0 && 2 != e1)
    return 2;
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "illegal"
        );
}
// =============================================================================
std::vector<int>
MeshTri::
computeBoundaryNodes_() const
{
  //TEUCHOS_ASSERT(false);
  return std::vector<int>();
}
// =============================================================================
}  // namespace Nosh
