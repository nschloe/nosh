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

#ifndef VIO_EPETRAMESH_MESH_H
#define VIO_EPETRAMESH_MESH_H

#include "VIO_Typedefs.h"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_SerialDenseVector.hpp>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>

namespace VIO {
namespace EpetraMesh {
class Mesh
{
  public:

    // constructor
    Mesh( const Teuchos::RCP<const Epetra_Comm> & comm,
          const Teuchos::RCP<const Epetra_Map>  & nodesMap,
          const Teuchos::RCP<const Epetra_Map>  & complexValuesMap
        );

    // destructor
    ~Mesh();

    unsigned int
    getNumNodes() const;

    const Teuchos::RCP<const Epetra_Map>
    getNodesMap() const;

    const Teuchos::RCP<const Epetra_Map>
    getComplexValuesMap() const;

    const Teuchos::ArrayRCP<Point>
    getNodes() const;

    Teuchos::ArrayRCP<Point>
    getNodesNonConst();

    void
    setNodes( const Teuchos::ArrayRCP<Point> nodes );

    void
    setBoundaryNodes( const Teuchos::ArrayRCP<bool> );

    const Teuchos::ArrayRCP<const bool>
    getBoundaryNodes() const;

    const Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> >
    getElems() const;

    Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> >
    getElemsNonConst();

    void
    setElems( const Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> > elems );

    // these are (approximately?) the valid VTK element types
    enum ElementType { EDGE2,    // line
                       EDGE3,    // quadratic edge
                       TRI3,     // triangle
                       TRI6,     // quadratic triangle
                       QUAD4,    // quad
                       QUAD8,    // quadratic quad
                       QUAD9,    // biquadratid quad
                       TET4,     // tetra
                       TET10,    // quadratic tetra
                       HEX8,     // hexahedron
                       HEX20,    // quadratic hexahedron
                       PRISM6,   // wegde
                       PRISM15,  // higher order wedge
                       PYRAMID5, // pyramid
                       INVALID
                     };

    void
    setElemTypes( const Teuchos::ArrayRCP<Mesh::ElementType> elemTypes );

    const Teuchos::ArrayRCP<const Mesh::ElementType>
    getElemTypes() const;

    Teuchos::RCP<Epetra_Vector>
    getControlVolumes() const;

    Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
    getEdgeLengths() const;

    Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
    getCoedgeLengths() const;

    double
    getDomainArea() const;

    void
    scale( double alpha );

    void
    scale( const Teuchos::Tuple<double,3> & alpha );

  protected:

  private:

      void
      computeFvmEntities_() const;

      Teuchos::RCP<Epetra_Map>
      getElemsToNodesMap_() const;

      double
      computeDomainArea_() const;

      Point
      add_( double alpha, const Point & x,
            double beta,  const Point & y
          ) const;

      //! Compute the area of a triangle given by
      //! \c x0, \c x1, and \c x2.
      double
      getTriangleArea_( const Point & x0,
                        const Point & x1,
                        const Point & x2
                      ) const;

      //! Compute the circumcenter of a triangle given by
      //! \c x0, \c x1, and \c x2.
      Point
      computeCircumcenter_( const Point & x0, const Point & x1, const Point & x2
                          ) const;

      double
      dot_( const Point & v, const Point & w
          ) const;

      //! Compute the cross product of two vectors
      //! \c v and \c w.
      Point
      cross_( const Point & v, const Point & w
            ) const;

      double
      norm2_( const Point & x
            ) const;

      //! Takes a vector \c x with an overlay map and makes sure that
      //! all the overlapping entries are summed up across all processes.
      void
      sumInOverlapMap_( Teuchos::RCP<Epetra_Vector> x ) const;

  private:

    const Teuchos::RCP<const Epetra_Comm> & comm_;
    const Teuchos::RCP<const Epetra_Map>  & nodesMap_;
    const Teuchos::RCP<const Epetra_Map>  & complexValuesMap_;

    Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> > elems_;
    Teuchos::ArrayRCP<ElementType> elemTypes_;
    Teuchos::ArrayRCP<Point> nodes_;
    Teuchos::ArrayRCP<bool> isBoundaryNode_;
    mutable Teuchos::RCP<Epetra_Vector> controlVolumes_;
    mutable Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeLengths_;
    mutable Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > coedgeLengths_;
    mutable double area_;
    Teuchos::Tuple<double,3> scaling_;
    mutable bool fvmEntitiesUpToDate_;
};
}
}

#endif // VIO_EPETRAMESH_MESH_H
