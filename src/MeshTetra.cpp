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
#include "MeshTetra.hpp"

#include <memory>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace Nosh
{
// =============================================================================
MeshTetra::
MeshTetra(
    const std::shared_ptr<const Teuchos::Comm<int>> & _comm,
    const std::shared_ptr<stk::io::StkMeshIoBroker> & broker
    ) :
  Mesh(_comm, broker),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  computeEdgeCoefficientsTime_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: MeshTetra::computeEdgeCoefficients"
        )),
  computeControlVolumesTime_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: MeshTetra::computeControlVolumes"
        )),
  computeBoundaryNodesTime_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: MeshTetra::computeBoundaryNodes"
        )),
#endif
  controlVolumes_(this->computeControlVolumes_()),
  edgeCoefficients_(this->computeEdgeCoefficients_()),
  boundaryNodeGids_(this->computeBoundaryNodeGids_())
{
}
// =============================================================================
MeshTetra::
~MeshTetra()
{
}
// =============================================================================
std::vector<double>
MeshTetra::
computeEdgeCoefficients_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
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

    auto edgeCoeffs = getEdgeCoefficientsCell_(localEdgeCoords);

    // Fill the edge coefficients into the vector.
    for (size_t i = 0; i < numLocalEdges; i++) {
      edgeCoefficients[edgeData_.cellEdges[k][i]] += edgeCoeffs[i];
    }
  }

  return edgeCoefficients;
}
// =============================================================================
Eigen::VectorXd
MeshTetra::
getEdgeCoefficientsCell_(
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
  double vol;
  // TODO Come up with a cleaner solution here.
  try {
    vol = getTetrahedronVolume_(edges[0], edges[1], edges[2]);
  } catch(...) {
    // If computing the volume throws an exception, then the edges chosen
    // happened to be conplanar. Changing one of those fixes this.
    vol = getTetrahedronVolume_(edges[0], edges[1], edges[3]);
  }

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
MeshTetra::
computeControlVolumes_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*computeControlVolumesTime_);
#endif
#ifndef NDEBUG
  TEUCHOS_ASSERT(nodesMap_);
  TEUCHOS_ASSERT(nodesOverlapMap_);
  TEUCHOS_ASSERT(ioBroker_);
#endif

  auto controlVolumes = std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(nodesMap_)
      );

  // Create temporaries to hold the overlap values for control volumes and
  // average thickness.
  Tpetra::Vector<double,int,int> cvOverlap(Teuchos::rcp(nodesOverlapMap_));

  this->computeControlVolumesT_(cvOverlap);

  // Export control volumes to a non-overlapping map, and sum the entries.
  Teuchos::RCP<const Tpetra::Export<int,int>> exporter = Tpetra::createExport(
      Teuchos::rcp(nodesOverlapMap_),
      Teuchos::rcp(nodesMap_)
      );
  controlVolumes->doExport(cvOverlap, *exporter, Tpetra::ADD);

  return controlVolumes;
}
// =============================================================================
void
MeshTetra::
computeControlVolumesT_(Tpetra::Vector<double,int,int> & cvOverlap) const
{
  std::vector<stk::mesh::Entity> cells = this->getOwnedCells();
  const size_t numCells = cells.size();

  const VectorFieldType & coordsField = getNodeField("coordinates");

  Teuchos::ArrayRCP<double> cvData = cvOverlap.getDataNonConst();

  // Calculate the contributions to the finite volumes cell by cell.
  for (size_t k = 0; k < numCells; k++) {
    const stk::mesh::Entity * localNodes =
      ioBroker_->bulk_data().begin_nodes(cells[k]);
    unsigned int numLocalNodes = ioBroker_->bulk_data().num_nodes(cells[k]);
#ifndef NDEBUG
    // Confirm that we always have the same simplices.
    TEUCHOS_ASSERT_EQUALITY(numLocalNodes, 4);
#endif

    // Fetch the nodal positions into 'localNodes'.
    std::vector<Eigen::Vector3d> localNodeCoords(numLocalNodes);
    for (unsigned int i = 0; i < numLocalNodes; i++) {
      localNodeCoords[i] = this->getNodeValue(coordsField, localNodes[i]);
    }

    // compute the circumcenter of the cell
    const auto cc = this->computeTetrahedronCircumcenter_(localNodeCoords);

    // Iterate over the edges.
    // As true edge entities are not available here, loop over all pairs of
    // local nodes.
    for (unsigned int e0 = 0; e0 < numLocalNodes; e0++) {
      const Eigen::Vector3d &x0 = localNodeCoords[e0];
      // TODO check if "- 1" is still needed
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

        // Get the other nodes.
        std::set<unsigned int> otherSet = this->getOtherIndices_(e0, e1);
        // Convert to vector (easier to handle for now)
        std::vector<unsigned int> other(otherSet.begin(), otherSet.end() );

        double edgeLength = (x1 - x0).norm();

        // Compute the (n-1)-dimensional covolume.
        const Eigen::Vector3d &other0 = localNodeCoords[other[0]];
        const Eigen::Vector3d &other1 = localNodeCoords[other[1]];
        double covolume = this->computeCovolume3d_(cc, x0, x1, other0, other1);
        // Throw an exception for 3D volumes.
        // To compute the average of the thicknesses of a control volume, one
        // has to loop over all the edges and add the thickness value to both
        // endpoints.  Then eventually, for each node, divide the resulting sum
        // by the number of connections (=number of faces of the finite
        // volume).  However, looping over edges is not (yet) possible. Hence,
        // we loop over all the cells here. This way, the edges are counted
        // several times, but it is difficult to determine how many times
        // exactly.
        //TEUCHOS_TEST_FOR_EXCEPTION(
        //    true,
        //    std::runtime_error,
        //    "Cannot calculate the average thickness in a 3D control volume yet."
        //    );

        // Compute the contributions to the finite volumes of the adjacent
        // edges.
        double pyramidVolume = 0.5*edgeLength * covolume / 3;
        cvData[lid0] += pyramidVolume;
        cvData[lid1] += pyramidVolume;
      }
    }
  }
}
// =============================================================================
double
MeshTetra::
computeCovolume3d_(
    const Eigen::Vector3d &cc,
    const Eigen::Vector3d &x0,
    const Eigen::Vector3d &x1,
    const Eigen::Vector3d &other0,
    const Eigen::Vector3d &other1
    ) const
{
  double covolume = 0.0;

  // edge midpoint
  Eigen::Vector3d mp = 0.5 * (x0 + x1);

  // Compute the circumcenters of the adjacent faces.
  // This could be precomputed as well.
  Eigen::Vector3d ccFace0 = computeTriangleCircumcenter_(x0, x1, other0);
  Eigen::Vector3d ccFace1 = computeTriangleCircumcenter_(x0, x1, other1);

  // Compute the area of the quadrilateral.
  // There are some really tricky degenerate cases here, i.e., combinations
  // of when ccFace{0,1}, cc, sit outside of the tetrahedron.

  // Use the triangle (MP, localNodes[other[0]], localNodes[other[1]]) (in this
  // order) to gauge the orientation of the two triangles that compose the
  // quadrilateral.
  Eigen::Vector3d gauge = (other0 - mp).cross(other1 - mp);

  // Add the area of the first triangle (MP,ccFace0,cc).
  // This makes use of the right angles.
  double triangleHeight0 = (mp - ccFace0).norm();
  double triangleArea0 = 0.5 * triangleHeight0 * (ccFace0 - cc).norm();

  // Check if the orientation of the triangle (MP,ccFace0,cc) coincides with
  // the orientation of the gauge triangle. If yes, add the area, subtract
  // otherwise.
  Eigen::Vector3d triangleNormal0 = (ccFace0 - mp).cross(cc - mp);

  // copysign takes the absolute value of the first argument and the sign of
  // the second.
  covolume += copysign(triangleArea0, triangleNormal0.dot(gauge));

  // Add the area of the second triangle (MP,cc,ccFace1).
  // This makes use of the right angles.
  double triangleHeight1 = (mp - ccFace1).norm();
  double triangleArea1 = 0.5 * triangleHeight1 * (ccFace1 - cc).norm();

  // Check if the orientation of the triangle (MP,cc,ccFace1) coincides with
  // the orientation of the gauge triangle. If yes, add the area, subtract
  // otherwise.
  Eigen::Vector3d triangleNormal1 = (cc - mp).cross(ccFace1 - mp);

  // copysign takes the absolute value of the first argument and the sign of
  // the second.
  covolume += copysign(triangleArea1, triangleNormal1.dot(gauge));

  return covolume;
}
// =============================================================================
std::set<unsigned int>
MeshTetra::
getOtherIndices_(unsigned int e0, unsigned int e1) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_INEQUALITY(e0, !=, e1);
#endif
  // Get the two indices in [0,1,2,3] which are not e0, e1.
  int myint[] = {0, 1, 2, 3};
  std::set<unsigned int> a(myint, myint+4);
  a.erase(e0);
  a.erase(e1);
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(a.size(), 2);
#endif
  return a;
}
// =============================================================================
double
MeshTetra::
getTetrahedronVolume_(
    const Eigen::Vector3d &edge0,
    const Eigen::Vector3d &edge1,
    const Eigen::Vector3d &edge2
    ) const
{
  // Make sure the edges are not conplanar.
  double alpha = edge0.dot(edge1.cross(edge2));
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      fabs(alpha) / edge0.norm() / edge1.norm() / edge2.norm() < 1.0e-5,
      "Illegal mesh: tetrahedron too flat.\n"
      << "The following edges (with origin (0,0,0)) "
      << "seem to be conplanar:\n\n"
      << "  (0) " << edge0 << ",\n"
      << "  (1) " << edge1 << ",\n"
      << "  (2) " << edge2 << ",\n\n"
      << "Abort."
      );
  double vol = fabs(alpha) / 6.0;
  return vol;
}
// =============================================================================
Eigen::Vector3d
MeshTetra::
computeTetrahedronCircumcenter_(
  const std::vector<Eigen::Vector3d> &nodes
  ) const
{
  // http://www.cgafaq.info/wiki/Tetrahedron_Circumsphere
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(nodes.size(), 4);
#endif

  // Compute with respect to the first point.
  std::vector<Eigen::Vector3d> relNodes(3);
  for (int k = 0; k < 3; k++) {
    relNodes[k] = nodes[k+1] - nodes[0];
  }

  double omega = 2.0 * relNodes[0].dot(relNodes[1].cross(relNodes[2]));

  // don't divide by 0
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      fabs(omega) < 1.0e-10,
      "It seems that the nodes \n"
      << "\n"
      << "   " << nodes[0] << "\n"
      << "   " << nodes[1] << "\n"
      << "   " << nodes[2] << "\n"
      << "   " << nodes[3] << "\n"
      << "\n"
      << "do not form a proper tetrahedron. Abort."
      << std::endl
      );
  const double alpha = relNodes[0].squaredNorm() / omega;
  const double beta  = relNodes[1].squaredNorm() / omega;
  const double gamma = relNodes[2].squaredNorm() / omega;

  return nodes[0]
    + alpha * relNodes[1].cross(relNodes[2])
    + beta *  relNodes[2].cross(relNodes[0])
    + gamma * relNodes[0].cross(relNodes[1]);
}
// =============================================================================
std::vector<stk::mesh::Entity>
MeshTetra::
getOverlapFaces_() const
{
  stk::mesh::Selector select_overlap_in_part =
    stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().universal_part())
    & (stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().locally_owned_part())
       |stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().globally_shared_part()));

  std::vector<stk::mesh::Entity> faces;
  stk::mesh::get_selected_entities(
      select_overlap_in_part,
      ioBroker_->bulk_data().buckets(stk::topology::FACE_RANK),
      faces
      );
  return faces;
}
// =============================================================================
std::set<int>
MeshTetra::
computeBoundaryNodeGids_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*computeBoundaryNodesTime_);
#endif

  auto myFaces = this->getOverlapFaces_();

  std::set<int> boundaryNodeGids;
  for (size_t k = 0; k < myFaces.size(); k++) {
    // if the face has one element, it's on the boundary
    if (ioBroker_->bulk_data().num_elements(myFaces[k]) == 1) {
      stk::mesh::Entity const * nodes =
        ioBroker_->bulk_data().begin_nodes(myFaces[k]);
#ifndef NDEBUG
      TEUCHOS_ASSERT_EQUALITY(ioBroker_->bulk_data().num_nodes(myFaces[k]), 3);
#endif
      boundaryNodeGids.insert(this->gid(nodes[0]));
      boundaryNodeGids.insert(this->gid(nodes[1]));
      boundaryNodeGids.insert(this->gid(nodes[2]));
    }
  }

  return boundaryNodeGids;
}
// =============================================================================
}  // namespace Nosh
