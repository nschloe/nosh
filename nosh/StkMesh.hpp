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

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Teuchos_Array.hpp>
#include <Teuchos_SerialDenseVector.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>
// =============================================================================
// forward declarations
namespace stk_classic
{
namespace mesh
{
class MetaData;
class BulkData;
}
namespace io
{
class MeshData;
}
} // namespace stk
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Map;
// =============================================================================
// typedefs
typedef stk_classic::mesh::Field<double, stk_classic::mesh::Cartesian> VectorFieldType;
typedef stk_classic::mesh::Field<double>                      ScalarFieldType;
typedef stk_classic::mesh::Field<int>                         IntScalarFieldType;
typedef Teuchos::SerialDenseVector<int, double>        DoubleVector;
typedef Teuchos::SerialDenseVector<int, const double>  ConstDoubleVector;
// =============================================================================
namespace Nosh
{

class StkMesh
{
private:
// Keep bulkData a pointer since its copy constructor is private; this causes
// issues when trying to copy (or initialize) MeshDataContainer.
  struct MeshDataContainer {
    stk_classic::mesh::fem::FEMMetaData metaData;
    Teuchos::RCP<stk_classic::io::MeshData> meshData;
    Teuchos::RCP<stk_classic::mesh::BulkData> bulkData;
  };

  struct EdgesContainer {
    //! Local edge ID -> Global node IDs.
    Teuchos::Array<Teuchos::Tuple<stk_classic::mesh::Entity*, 2> > edgeNodes;
    //! Local cell ID -> Local edge IDs.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> > cellEdges;
  };

public:
  StkMesh(const Epetra_Comm &comm,
          const std::string &fileName,
          const int index
        );

  virtual
  ~StkMesh();

  MeshDataContainer
  read(const Epetra_Comm &comm,
       const std::string &fileName,
       const int index
     );

  double
  getTime() const;

  void
  openOutputChannel(const std::string &outputDir,
                    const std::string &fileBaseName
                   );

  void
  write(const Epetra_Vector & psi,
        const double time
       ) const;

  Teuchos::RCP<Epetra_Vector>
  createVector(const std::string & fieldName) const;

  Teuchos::RCP<Epetra_MultiVector>
  createMultiVector(const std::string & fieldName) const;

  Teuchos::RCP<Epetra_Vector>
  createComplexVector(const std::string & fieldName) const;

  void
  mergeComplexVector_(const Epetra_Vector &psi,
                      const std::string & fieldName
                    ) const;

  unsigned int
  getNumNodes() const;

  Teuchos::RCP<const Epetra_Vector>
  getControlVolumes() const;

  double
  getDomainVolume() const;

  const Epetra_Comm &
  getComm() const;

  Teuchos::ArrayRCP<double>
  getEdgeCoefficients() const;

  std::vector<stk_classic::mesh::Entity*>
  getOwnedCells() const;

  std::vector<stk_classic::mesh::Entity*>
  getOverlapEdges() const;

  const Teuchos::Array<Teuchos::Tuple<stk_classic::mesh::Entity*, 2> >
  getEdgeNodes() const;

  std::vector<stk_classic::mesh::Entity*>
  getOwnedNodes() const;

  std::vector<stk_classic::mesh::Entity*>
  getOverlapNodes() const;

//const DoubleVector
//getNodeCoordinatesNonconst(const stk_classic::mesh::Entity * nodeEntity) const;

  Teuchos::RCP<const Epetra_Map>
  getNodesMap() const;

  Teuchos::RCP<const Epetra_Map>
  getNodesOverlapMap() const;

  Teuchos::RCP<const Epetra_Map>
  getComplexNonOverlapMap() const;

  Teuchos::RCP<const Epetra_Map>
  getComplexOverlapMap() const;

  const Epetra_Map&
  getComplexNonOverlapMap2() const;

  unsigned int
  getNumEdgesPerCell(unsigned int cellDimension) const;

  const DoubleVector
  getVectorFieldNonconst(const stk_classic::mesh::Entity * nodeEntity,
                         const std::string & fieldName,
                         const int numDims
                        ) const;
  double
  getScalarFieldNonconst(const stk_classic::mesh::Entity * nodeEntity,
                         const std::string & fieldName
                        ) const;

protected:
private:
  const int numDim_;

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> computeEdgeCoefficientsTime_;
  const Teuchos::RCP<Teuchos::Time> writeTime_;
#endif

  const Epetra_Comm &comm_;

  MeshDataContainer meshDataContainer_;

  const std::vector<stk_classic::mesh::Entity*> ownedNodes_;

  const Teuchos::RCP<const Epetra_Map> nodesMap_;
  const Teuchos::RCP<const Epetra_Map> nodesOverlapMap_;
  const Teuchos::RCP<const Epetra_Map> complexMap_;
  const Teuchos::RCP<const Epetra_Map> complexOverlapMap_;

  const Teuchos::RCP<const Epetra_Vector> controlVolumes_;

  const EdgesContainer edgeData_;

  const Teuchos::ArrayRCP<double> edgeCoefficients_;

  bool outputChannelIsOpen_;

  double time_;

private:
  void
  computeControlVolumesTri_(
      const Teuchos::RCP<Epetra_Vector> & cvOverlap
      ) const;

  void
  computeControlVolumesTet_(
      const Teuchos::RCP<Epetra_Vector> & cvOverlap
      ) const;

  Teuchos::RCP<Epetra_Vector>
  complexfield2vector_(const ScalarFieldType &realField,
                       const ScalarFieldType &imagField
                     ) const;

  Teuchos::RCP<Epetra_Vector>
  field2vector_(const ScalarFieldType &field) const;

  Teuchos::RCP<Epetra_MultiVector>
  field2vector_(const VectorFieldType &field,
                const int numComponents
               ) const;

  std::vector<stk_classic::mesh::Entity*>
  buildOwnedNodes_() const;

//Teuchos::ArrayRCP<const DoubleVector>
//getNodeCoordinates_(const stk_classic::mesh::PairIterRelation &relation) const;

  Teuchos::ArrayRCP<double>
  computeEdgeCoefficients_() const;

//! Compute the volume of the (Voronoi) control cells for each point.
  Teuchos::RCP<Epetra_Vector>
  computeControlVolumes_() const;

  Teuchos::RCP<const Epetra_Map>
  createEntitiesMap_(const std::vector<stk_classic::mesh::Entity*> &entityList) const;

  Teuchos::RCP<const Epetra_Map>
  createComplexMap_(const std::vector<stk_classic::mesh::Entity*> &nodeList) const;

  double
  computeCovolume2d_(const DoubleVector &cc,
                     const DoubleVector &x0,
                     const DoubleVector &x1,
                     const DoubleVector &other0
                    ) const;

  double
  computeCovolume3d_(const DoubleVector &cc,
                     const DoubleVector &x0,
                     const DoubleVector &x1,
                     const DoubleVector &other0,
                     const DoubleVector &other1
                    ) const;

  Teuchos::Tuple<unsigned int, 2>
  getOtherIndices_(unsigned int e0, unsigned int e1) const;

  DoubleVector
  add_(double alpha, const DoubleVector &x,
       double beta,  const DoubleVector &y
      ) const;


  double
  getTriangleArea_(const DoubleVector &edge0,
                   const DoubleVector &edge1
                  ) const;

  double
  getTetrahedronVolume_(const DoubleVector &edge0,
                        const DoubleVector &edge1,
                        const DoubleVector &edge2
                       ) const;

  DoubleVector
  computeTriangleCircumcenter_(const DoubleVector &node0,
                               const DoubleVector &node1,
                               const DoubleVector &node2
                              ) const;
  DoubleVector
  computeTriangleCircumcenter_(
    const Teuchos::ArrayRCP<const DoubleVector> &nodes
    ) const;

  DoubleVector
  computeTetrahedronCircumcenter_(
    const Teuchos::ArrayRCP<const DoubleVector> &nodes
    ) const;

  DoubleVector
  getEdgeCoefficientsNumerically_(
    const Teuchos::ArrayRCP<const DoubleVector> edges
    ) const;

  double
  dot_(const DoubleVector &v, const DoubleVector &w) const;

  DoubleVector
  cross_(const DoubleVector &v,
         const DoubleVector &w
        ) const;

  double
  norm2_(const DoubleVector &x) const;

  double
  norm2squared_(const DoubleVector &x) const;

  EdgesContainer
  createEdgeData_();
};
// -----------------------------------------------------------------------------
} // namespace Nosh
// =============================================================================
#endif // NOSH_STKMESH_H
