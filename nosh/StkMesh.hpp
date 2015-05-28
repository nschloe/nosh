// @HEADER
//
//    Mesh class with compatibility to stk_mesh.
//    Copyright (C) 2010--2012  Nico Schl\"omer
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
#ifndef NOSH_STKMESH_H
#define NOSH_STKMESH_H
// =============================================================================
// includes
#include <string>
#include <vector>
#include <tuple>
#include <set>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Epetra_IntSerialDenseVector.h>
#include <Epetra_FECrsGraph.h>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <Eigen/Dense>
// =============================================================================
// forward declarations
namespace stk
{
namespace mesh
{
class BulkData;
}
} // namespace stk
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Map;
// =============================================================================
// typedefs
typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::mesh::Field<int> IntScalarFieldType;
typedef std::tuple<stk::mesh::Entity, stk::mesh::Entity> edge;
// =============================================================================
namespace Nosh
{

class StkMesh
{
private:
  // Keep bulkData a pointer since its copy constructor is private; this causes
  // issues when trying to copy (or initialize) MeshDataContainer.
  struct EdgesContainer {
    //! Local edge ID -> Global node IDs.
    std::vector<edge> edgeNodes;
    //! Local cell ID -> Local edge IDs.
    std::vector<std::vector<int> > cellEdges;
  };

public:
  StkMesh(
      const std::shared_ptr<const Epetra_Comm> & comm,
      const std::string &fileName,
      const int index
      );

  virtual
  ~StkMesh();

  double
  getTime() const;

  void
  openOutputChannel(const std::string &outputDir,
                    const std::string &fileBaseName
                   );

  void
  write(const double time) const;

  std::shared_ptr<Epetra_Vector>
  createVector(const std::string & fieldName) const;

  std::shared_ptr<Epetra_MultiVector>
  createMultiVector(const std::string & fieldName) const;

  std::shared_ptr<Epetra_Vector>
  createComplexVector(const std::string & fieldName) const;

  void
  insert(const Epetra_Vector & psi,
         const std::string & fieldName
         ) const;

  unsigned int
  getNumNodes() const;

  std::shared_ptr<const Epetra_Vector>
  getControlVolumes() const;

  double
  getDomainVolume() const;

  const Epetra_Comm &
  getComm() const;

  std::vector<double>
  getEdgeCoefficients() const;

  std::vector<stk::mesh::Entity>
  getOwnedCells() const;

  std::vector<stk::mesh::Entity>
  getOverlapEdges() const;

  const std::vector<std::tuple<stk::mesh::Entity, stk::mesh::Entity> >
  getEdgeNodes() const;

  std::vector<stk::mesh::Entity>
  getOwnedNodes() const;

  std::vector<stk::mesh::Entity>
  getOverlapNodes() const;

//const Eigen::Vector3d
//getNodeCoordinatesNonconst(stk::mesh::Entity nodeEntity) const;

  std::shared_ptr<const Epetra_Map>
  getNodesMap() const;

  std::shared_ptr<const Epetra_Map>
  getNodesOverlapMap() const;

  std::shared_ptr<const Epetra_Map>
  getComplexNonOverlapMap() const;

  std::shared_ptr<const Epetra_Map>
  getComplexOverlapMap() const;

  const Epetra_Map&
  getComplexNonOverlapMap2() const;

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
  gid(const stk::mesh::Entity e) const;

  const Epetra_FECrsGraph
  buildComplexGraph() const;

protected:
private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const std::shared_ptr<Teuchos::Time> computeEdgeCoefficientsTime_;
  const std::shared_ptr<Teuchos::Time> writeTime_;
#endif

  const std::shared_ptr<const Epetra_Comm> comm_;

  // Apparently, process_output_request is not const. Make
  // the ioBroker_ mutable so our write() can be const.
  const std::shared_ptr<stk::io::StkMeshIoBroker> ioBroker_;

  const std::vector<stk::mesh::Entity> ownedNodes_;

  const std::shared_ptr<const Epetra_Map> nodesMap_;
  const std::shared_ptr<const Epetra_Map> nodesOverlapMap_;
  const std::shared_ptr<const Epetra_Map> complexMap_;
  const std::shared_ptr<const Epetra_Map> complexOverlapMap_;

  const std::shared_ptr<const Epetra_Vector> controlVolumes_;

  const EdgesContainer edgeData_;

  const std::vector<double> edgeCoefficients_;

  int outputChannel_;

  double time_;

public:
  const std::vector<Epetra_IntSerialDenseVector> globalIndexCache;

private:

  const std::vector<Epetra_IntSerialDenseVector>
  buildGlobalIndexCache_() const;

  std::shared_ptr<stk::io::StkMeshIoBroker>
  read_(const std::string &fileName,
        const int index
        );

  void
  computeControlVolumesTri_(Epetra_Vector & cvOverlap) const;

  void
  computeControlVolumesTet_(Epetra_Vector & cvOverlap) const;

  std::shared_ptr<Epetra_Vector>
  complexfield2vector_(const ScalarFieldType &realField,
                       const ScalarFieldType &imagField
                       ) const;

  std::shared_ptr<Epetra_Vector>
  field2vector_(const ScalarFieldType &field) const;

  std::shared_ptr<Epetra_MultiVector>
  field2vector_(const VectorFieldType &field,
                const int numComponents
                ) const;

  void
  mergeComplexVector_(const Epetra_Vector &psi,
                      const std::string & fieldName
                      ) const;

  std::vector<stk::mesh::Entity>
  buildOwnedNodes_(const stk::mesh::BulkData & myBulkData) const;

  std::vector<double>
  computeEdgeCoefficients_() const;

  //! Compute the volume of the (Voronoi) control cells for each point.
  std::shared_ptr<Epetra_Vector>
  computeControlVolumes_() const;

  std::shared_ptr<const Epetra_Map>
  createEntitiesMap_(const std::vector<stk::mesh::Entity> &entityList) const;

  std::shared_ptr<const Epetra_Map>
  createComplexMap_(const std::vector<stk::mesh::Entity> &nodeList) const;

  double
  computeCovolume2d_(const Eigen::Vector3d &cc,
                     const Eigen::Vector3d &x0,
                     const Eigen::Vector3d &x1,
                     const Eigen::Vector3d &other0
                    ) const;

  double
  computeCovolume3d_(const Eigen::Vector3d &cc,
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
StkMesh
read(const std::string & fileName,
     const int index
     );
} // namespace Nosh
// =============================================================================
#endif // NOSH_STKMESH_H
