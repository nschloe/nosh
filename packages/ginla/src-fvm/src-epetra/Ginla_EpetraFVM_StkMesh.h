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
#include "Ginla_Typedefs.h"

#include <string>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
// =============================================================================
// forward declarations
namespace stk {
    namespace mesh {
        class MetaData;
        class BulkData;
    }
}
class Epetra_Vector;
class Epetra_Map;
// =============================================================================
// typedefs
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;
// =============================================================================
namespace Ginla {
namespace EpetraFVM {

class StkMesh
{
public:
    StkMesh( const Epetra_Comm                       & comm,
             const Teuchos::RCP<stk::mesh::MetaData> & metaData,
             const Teuchos::RCP<stk::mesh::BulkData> & bulkData,
             const Teuchos::RCP<VectorFieldType>     & coordinatesField
           );

    virtual
    ~StkMesh();

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
    getCoedgeEdgeRatios() const;

    std::vector<stk::mesh::Entity*>
    getOwnedCells() const;

    Teuchos::Tuple<Point,3>
    getNodeCoordinates( const stk::mesh::PairIterRelation & relation ) const;

    Teuchos::RCP<Epetra_Map>
    getComplexMap() const;

protected:
private:

    const Epetra_Comm & comm_;

    const Teuchos::RCP<stk::mesh::MetaData> metaData_;
    const Teuchos::RCP<stk::mesh::BulkData> bulkData_;
    const Teuchos::RCP<VectorFieldType>     coordinatesField_;

    const Teuchos::RCP<Epetra_Map> nodesMap_;
    const Teuchos::RCP<Epetra_Map> nodesOverlapMap_;
    const Teuchos::RCP<Epetra_Map> complexMap_;

    Teuchos::Tuple<double,3> scaling_;

    mutable bool fvmEntitiesUpToDate_;

    const Teuchos::RCP<Epetra_Vector> controlVolumes_;
    const Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > coedgeEdgeRatio_;

    mutable double area_;

private:
    std::vector<stk::mesh::Entity*>
    getOwnedNodes_() const;

    std::vector<stk::mesh::Entity*>
    getOverlapNodes_() const;

    Teuchos::RCP<Epetra_Map>
    createMap_( const std::vector<stk::mesh::Entity*> & nodesList ) const;

    Teuchos::RCP<Epetra_Map>
    createComplexMap_() const;

    double
    computeDomainArea_() const;

    void
    computeFvmEntities_() const;

    Point
    add_( double alpha, const Point & x,
          double beta,  const Point & y
        ) const;

    double
    getTriangleArea_( const Point & x0,
                      const Point & x1,
                      const Point & x2
                    ) const;

    Point
    computeCircumcenter_( const Teuchos::Tuple<Point,3> & nodes ) const;

    double
    dot_( const Point & v, const Point & w
        ) const;

    Point
    cross_( const Point & v, const Point & w
          ) const;

    double
    norm2_( const Point & x
          ) const;
};
} // namespace EpetraFVM
} // namespace Ginla
// =============================================================================
#endif // GINLA_EPETRAFVM_STKMESH_H
