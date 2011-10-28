/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef GINLA_EPETRAFVM_STKMESH_H
#define GINLA_EPETRAFVM_STKMESH_H
// =============================================================================
// includes
// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <string>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>
// =============================================================================
// forward declarations
namespace stk {
    namespace mesh {
        class MetaData;
        class BulkData;
    }
    namespace io {
        class MeshData;
    }
}
class Epetra_Vector;
class Epetra_Map;
// =============================================================================
// typedefs
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;
typedef Teuchos::Tuple<double,3>                      Point;
// =============================================================================
namespace Ginla {
namespace EpetraFVM {

class StkMesh
{
public:
    StkMesh( const Epetra_Comm                               & comm,
             const Teuchos::RCP<stk::mesh::fem::FEMMetaData> & metaData,
             const Teuchos::RCP<stk::mesh::BulkData>         & bulkData,
             const Teuchos::RCP<VectorFieldType>             & coordinatesField
           );

    virtual
    ~StkMesh();

    void
    setOutputFile( const string & outputDir,
                   const string & fileBaseName
                 );


    const Teuchos::RCP<stk::mesh::fem::FEMMetaData>
    getMetaData() const;

    const Teuchos::RCP<stk::io::MeshData>
    getMeshData() const;

    const Teuchos::RCP<stk::mesh::BulkData>
    getBulkData() const;

    unsigned int
    getNumNodes() const;

    Teuchos::RCP<const Epetra_Vector>
    getControlVolumes() const;

    double
    getDomainArea() const;

    const Epetra_Comm &
    getComm() const;

    void
    scale( const Teuchos::Tuple<double,3> & newScaling );

    Teuchos::Tuple<double,3>
    getScaling() const;

    Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
    getEdgeCoefficients() const;

    std::vector<stk::mesh::Entity*>
    getOwnedCells() const;

    std::vector<stk::mesh::Entity*>
    getOwnedEdges() const;

    std::vector<stk::mesh::Entity*>
    getOwnedNodes() const;

    std::vector<stk::mesh::Entity*>
    getOverlapNodes() const;

    Teuchos::Array<Point>
    getNodeCoordinates( const stk::mesh::PairIterRelation & relation ) const;

    double
    getThickness( const stk::mesh::PairIterRelation & relation ) const;

    Teuchos::RCP<Epetra_Map>
    getNodesMap() const;

    Teuchos::RCP<Epetra_Map>
    getNodesOverlapMap() const;

    Teuchos::RCP<Epetra_Map>
    getComplexMap() const;

    Teuchos::RCP<Epetra_Map>
    getComplexOverlapMap() const;

    unsigned int
    getNumEdgesPerCell( unsigned int cellDimension ) const;

    unsigned int
    getCellDimension( const unsigned int numLocalNodes ) const;

    // Doesn't actually need to be public.
    void
    computeFvmEntities_() const;

protected:
private:

    const Epetra_Comm & comm_;

    const Teuchos::RCP<stk::mesh::fem::FEMMetaData> metaData_;
    const Teuchos::RCP<stk::io::MeshData>           meshData_;
    const Teuchos::RCP<stk::mesh::BulkData>         bulkData_;
    const Teuchos::RCP<VectorFieldType>             coordinatesField_;
//     const Teuchos::RCP<VectorFieldType>         thicknessField_;

    const Teuchos::RCP<Epetra_Map> nodesMap_;
    const Teuchos::RCP<Epetra_Map> nodesOverlapMap_;
    const Teuchos::RCP<Epetra_Map> complexMap_;
    const Teuchos::RCP<Epetra_Map> complexOverlapMap_;

    Teuchos::Tuple<double,3> scaling_;

    mutable bool fvmEntitiesUpToDate_;

    const Teuchos::RCP<Epetra_Vector> controlVolumes_;
    const Teuchos::RCP<Epetra_Vector> averageThickness_;
    const Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeCoefficients_;

    mutable double area_;

    bool isOutputFileSet_;

private:

    Teuchos::RCP<Epetra_Map>
    createNodesMap_( const std::vector<stk::mesh::Entity*> & nodesList ) const;

    Teuchos::RCP<Epetra_Map>
    createComplexMap_( const std::vector<stk::mesh::Entity*> & nodeList ) const;

    double
    computeCovolume2d_( const Point & cc,
                        const Point & x0,
                        const Point & x1,
                        const Point & other0
                      ) const;

    double
    computeCovolume3d_( const Point & cc,
                        const Point & x0,
                        const Point & x1,
                        const Point & other0,
                        const Point & other1
                      ) const;

    Teuchos::Tuple<unsigned int,2>
    getOtherIndices_( unsigned int e0, unsigned int e1 ) const;

    Point
    add_( double alpha, const Point & x,
          double beta,  const Point & y
        ) const;

    double
    getTriangleArea_( const Point & x0,
                      const Point & x1,
                      const Point & x2
                    ) const;

    double
    getTetrahedronVolume_( const Point & node0,
                           const Point & node1,
                           const Point & node2,
                           const Point & node3
                         ) const;

    Point
    computeTriangleCircumcenter_( const Point & node0, const Point & node1, const Point & node2 ) const;
    Point
    computeTriangleCircumcenter_( const Teuchos::Array<Point> & nodes ) const;

    Point
    computeTetrahedronCircumcenter_( const Teuchos::Array<Point> & nodes ) const;

    Teuchos::ArrayRCP<double>
    getEdgeCoefficientsNumerically_( const Teuchos::Array<Point> localNodes ) const;

    double
    dot_( const Point & v, const Point & w
        ) const;

    Point
    cross_( const Point & v, const Point & w
          ) const;

    double
    norm2_( const Point & x
          ) const;

    double
    norm2squared_( const Point & x
                ) const;
};
} // namespace EpetraFVM
} // namespace Ginla
// =============================================================================
#endif // GINLA_EPETRAFVM_STKMESH_H
