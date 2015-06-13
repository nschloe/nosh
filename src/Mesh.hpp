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
#ifndef NOSH_MESH_HPP
#define NOSH_MESH_HPP
// =============================================================================
// includes
#include <string>
#include <vector>
#include <tuple>
#include <set>

#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsGraph.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <Eigen/Dense>

// forward declarations
namespace stk
{
namespace mesh
{
class BulkData;
}
} // namespace stk

// typedefs
typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::mesh::Field<int> IntScalarFieldType;
typedef std::tuple<stk::mesh::Entity, stk::mesh::Entity> edge;

namespace Nosh
{

class Mesh
{
private:
  // Keep bulkData a pointer since its copy constructor is private; this causes
  // issues when trying to copy (or initialize) MeshDataContainer.
  struct EdgesContainer {
    //! Local edge ID -> Global node IDs.
    std::vector<edge> edgeNodes;
    //! Local cell ID -> Local edge IDs.
    std::vector<std::vector<int>> cellEdges;
  };

public:
  Mesh(
      const std::shared_ptr<const Teuchos::Comm<int>> & comm,
      const std::string &fileName,
      const int index
      );

  virtual
  ~Mesh();

  double
  getTime() const
  {
    return time_;
  };

  void
  openOutputChannel(const std::string &outputDir,
                    const std::string &fileBaseName
                   );

  void
  write(const double time) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  createVector(const std::string & fieldName) const;

  std::shared_ptr<Tpetra::MultiVector<double,int,int>>
  createMultiVector(const std::string & fieldName) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  createComplexVector(const std::string & fieldName) const;

  void
  insert(const Tpetra::Vector<double,int,int> & psi,
         const std::string & fieldName
         ) const;

  unsigned int
  getNumNodes() const
  {
    return nodesMap_->getGlobalNumElements();
  }

  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  getControlVolumes() const
  {
    return controlVolumes_;
  }

  double
  getDomainVolume() const
  {
    return controlVolumes_->norm1();
  }

  std::shared_ptr<const Teuchos::Comm<int>>
  getComm() const
  {
    return comm_;
  }

  std::vector<double>
  getEdgeCoefficients() const
  {
    return edgeCoefficients_;
  }

  std::vector<stk::mesh::Entity>
  getOwnedCells() const;

  std::vector<stk::mesh::Entity>
  getOverlapEdges() const;

  const std::vector<std::tuple<stk::mesh::Entity, stk::mesh::Entity>>
  getMyEdges() const
  {
    return edgeData_.edgeNodes;
  }

  std::vector<stk::mesh::Entity>
  getOwnedNodes() const
  {
    return ownedNodes_;
  }

  std::vector<stk::mesh::Entity>
  getOverlapNodes() const;

//const Eigen::Vector3d
//getNodeCoordinatesNonconst(stk::mesh::Entity nodeEntity) const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  getMap() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(nodesMap_);
#endif
    return nodesMap_;
  }

  std::shared_ptr<const Tpetra::Map<int,int>>
  getOverlapMap() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(nodesOverlapMap_);
#endif
    return nodesOverlapMap_;
  }

  std::shared_ptr<const Tpetra::Map<int,int>>
  getMapComplex() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(complexMap_);
#endif
    return complexMap_;
  }

  std::shared_ptr<const Tpetra::Map<int,int>>
  getOverlapMapComplex() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(complexOverlapMap_);
#endif
    return complexOverlapMap_;
  }

  unsigned int
  getNumEdgesPerCell(unsigned int cellDimension) const;

  const VectorFieldType &
  getNodeField(const std::string & fieldName) const;

  const Eigen::Vector3d
  getNodeValue(
      const VectorFieldType & field,
      stk::mesh::Entity nodeEntity
      ) const;

  double
  getScalarFieldNonconst(
      stk::mesh::Entity nodeEntity,
      const std::string & fieldName
      ) const;

  uint64_t
  gid(const stk::mesh::Entity e) const
  {
    return ioBroker_->bulk_data().identifier(e) - 1;
  }

  Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
  buildGraph() const;

  Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
  buildComplexGraph() const;

protected:
private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const std::shared_ptr<Teuchos::Time> computeEdgeCoefficientsTime_;
  const std::shared_ptr<Teuchos::Time> writeTime_;
#endif

  const std::shared_ptr<const Teuchos::Comm<int>> comm_;

  // Apparently, process_output_request is not const. Make
  // the ioBroker_ mutable so our write() can be const.
  const std::shared_ptr<stk::io::StkMeshIoBroker> ioBroker_;

  const std::vector<stk::mesh::Entity> ownedNodes_;

  const std::shared_ptr<const Tpetra::Map<int,int>> nodesMap_;
  const std::shared_ptr<const Tpetra::Map<int,int>> nodesOverlapMap_;
  const std::shared_ptr<const Tpetra::Map<int,int>> complexMap_;
  const std::shared_ptr<const Tpetra::Map<int,int>> complexOverlapMap_;

  const std::shared_ptr<const Tpetra::Vector<double,int,int>> controlVolumes_;

  const EdgesContainer edgeData_;

  const std::vector<double> edgeCoefficients_;

  int outputChannel_;

  double time_;

public:
  const std::vector<Teuchos::Tuple<int,2>> edgeGids;
  const std::vector<Teuchos::Tuple<int,4>> edgeGidsComplex;

private:

  const std::vector<Teuchos::Tuple<int,2>>
  buildEdgeGids_() const;

  const std::vector<Teuchos::Tuple<int,4>>
  buildEdgeGidsComplex_() const;

  std::shared_ptr<stk::io::StkMeshIoBroker>
  read_(
      const std::string &fileName,
      const int index
      );

  void
  computeControlVolumesTri_(Tpetra::Vector<double,int,int> & cvOverlap) const;

  void
  computeControlVolumesTet_(Tpetra::Vector<double,int,int> & cvOverlap) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  complexfield2vector_(
      const ScalarFieldType &realField,
      const ScalarFieldType &imagField
      ) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  field2vector_(const ScalarFieldType &field) const;

  std::shared_ptr<Tpetra::MultiVector<double,int,int>>
  field2vector_(
      const VectorFieldType &field,
      const int numComponents
      ) const;

  void
  mergeComplexVector_(
      const Tpetra::Vector<double,int,int> &psi,
      const std::string & fieldName
      ) const;

  std::vector<stk::mesh::Entity>
  buildOwnedNodes_(const stk::mesh::BulkData & myBulkData) const;

  std::vector<double>
  computeEdgeCoefficients_() const;

  //! Compute the volume of the (Voronoi) control cells for each point.
  std::shared_ptr<Tpetra::Vector<double,int,int>>
  computeControlVolumes_() const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  createEntitiesMap_(const std::vector<stk::mesh::Entity> &entityList) const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  createComplexMap_(const std::vector<stk::mesh::Entity> &nodeList) const;

  double
  computeCovolume2d_(
      const Eigen::Vector3d &cc,
      const Eigen::Vector3d &x0,
      const Eigen::Vector3d &x1,
      const Eigen::Vector3d &other0
      ) const;

  double
  computeCovolume3d_(
      const Eigen::Vector3d &cc,
      const Eigen::Vector3d &x0,
      const Eigen::Vector3d &x1,
      const Eigen::Vector3d &other0,
      const Eigen::Vector3d &other1
      ) const;

  unsigned int
  getOtherIndex_(unsigned int e0, unsigned int e1) const;

  std::set<unsigned int>
  getOtherIndices_(unsigned int e0, unsigned int e1) const;

  double
  getTriangleArea_(const Eigen::Vector3d &edge0,
                   const Eigen::Vector3d &edge1
                  ) const;

  double
  getTetrahedronVolume_(
      const Eigen::Vector3d &edge0,
      const Eigen::Vector3d &edge1,
      const Eigen::Vector3d &edge2
      ) const;

  Eigen::Vector3d
  computeTriangleCircumcenter_(
      const Eigen::Vector3d &node0,
      const Eigen::Vector3d &node1,
      const Eigen::Vector3d &node2
      ) const;

  Eigen::Vector3d
  computeTriangleCircumcenter_(
    const std::vector<Eigen::Vector3d> &nodes
    ) const;

  Eigen::Vector3d
  computeTetrahedronCircumcenter_(
    const std::vector<Eigen::Vector3d> &nodes
    ) const;

  Eigen::VectorXd
  getEdgeCoefficientsNumerically_(
    const std::vector<Eigen::Vector3d> edges
    ) const;

  double
  norm2squared_(const Eigen::Vector3d &x) const;

  EdgesContainer
  createEdgeData_();
};
// -----------------------------------------------------------------------------
// Helper function
Mesh
read(
    const std::string & fileName,
    const int index = 0
    );

} // namespace Nosh
// =============================================================================
#endif // NOSH_MESH_HPP
