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
#include <memory>
#include <string>
#include <vector>
#include <tuple>

#include <Teuchos_RCP.hpp>
//#ifdef NOSH_TEUCHOS_TIME_MONITOR
//#include <Teuchos_Time.hpp>
//#endif
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsGraph.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
//#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldTraits.hpp>

#include <Eigen/Dense>

// typedefs
typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;
//typedef stk::mesh::Field<int> IntScalarFieldType;
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
      const std::shared_ptr<stk::io::StkMeshIoBroker> & broker
      );

  virtual
  ~Mesh();

  double
  getTime() const
  {
    return time_;
  };

  void
  openFile(const std::string &outputFile);

  void
  write(const double time = 0.0) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  getVector(const std::string & fieldName) const;

  std::shared_ptr<Tpetra::MultiVector<double,int,int>>
  getMultiVector(const std::string & fieldName) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  getComplexVector(const std::string & fieldName) const;

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

  const VectorFieldType &
  getNodeField(const std::string & fieldName) const;

  const Eigen::Vector3d
  getNodeValue(
      const VectorFieldType & field,
      stk::mesh::Entity nodeEntity
      ) const
  {
    return Eigen::Vector3d(stk::mesh::field_data(field, nodeEntity));
  };

  uint64_t
  gid(const stk::mesh::Entity e) const
  {
    return ioBroker_->bulk_data().identifier(e) - 1;
  }

  Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
  buildGraph() const;

  Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
  buildComplexGraph() const;

  void
  insertVector(
      const Tpetra::Vector<double,int,int> &x,
      const std::string & fieldName
      ) const;

  void
  insertComplexVector(
      const Tpetra::Vector<double,int,int> &psi,
      const std::string & fieldName
      ) const;

public:
  virtual
  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  getControlVolumes() const = 0;

  virtual
  std::vector<double>
  getEdgeCoefficients() const = 0;

  virtual
  std::vector<int>
  getBoundaryNodes() const = 0;

protected:

  Eigen::Vector3d
  computeTriangleCircumcenter_(
      const std::vector<Eigen::Vector3d> &nodes
      ) const;

  Eigen::Vector3d
  computeTriangleCircumcenter_(
      const Eigen::Vector3d &node0,
      const Eigen::Vector3d &node1,
      const Eigen::Vector3d &node2
      ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const std::shared_ptr<Teuchos::Time> computeEdgeCoefficientsTime_;
  const std::shared_ptr<Teuchos::Time> writeTime_;
#endif

public:
  const std::shared_ptr<const Teuchos::Comm<int>> comm;

protected:
  // Apparently, process_output_request is not const. Make
  // the ioBroker_ mutable so our write() can be const.
  const std::shared_ptr<stk::io::StkMeshIoBroker> ioBroker_;

private:
  const std::vector<stk::mesh::Entity> ownedNodes_;

protected:
  const std::shared_ptr<const Tpetra::Map<int,int>> nodesMap_;
  const std::shared_ptr<const Tpetra::Map<int,int>> nodesOverlapMap_;

private:
  const std::shared_ptr<const Tpetra::Map<int,int>> complexMap_;
  const std::shared_ptr<const Tpetra::Map<int,int>> complexOverlapMap_;

protected:
  const EdgesContainer edgeData_;

private:
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

  std::vector<stk::mesh::Entity>
  buildOwnedNodes_(const stk::mesh::BulkData & myBulkData) const;

  std::vector<double>
  computeEdgeCoefficients_() const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  createEntitiesMap_(const std::vector<stk::mesh::Entity> &entityList) const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  createComplexMap_(const std::vector<stk::mesh::Entity> &nodeList) const;

  EdgesContainer
  createEdgeData_();
};
// -----------------------------------------------------------------------------

} // namespace Nosh
// =============================================================================
#endif // NOSH_MESH_HPP
