/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Sch\"omer

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

#include "VIO_Mesh_Mesh.h"
#include <Epetra_Map.h>

// =============================================================================
VIO::Mesh::Mesh::
Mesh( const Teuchos::RCP<const Teuchos::Comm<int> > & comm ):
    comm_( comm ),
    elems_( Teuchos::null ),
    elemTypes_( Teuchos::null ),
    nodes_( Teuchos::null ),
    isBoundaryNode_( Teuchos::null ),
    controlVolumes_( Teuchos::null ),
    area_( 0.0 )
{
}
// =============================================================================
VIO::Mesh::Mesh::
~Mesh()
{
}
// =============================================================================
unsigned int
VIO::Mesh::Mesh::
getNumNodes() const
{
  return nodes_.size();
}
// =============================================================================
const Teuchos::ArrayRCP<Point>
VIO::Mesh::Mesh::
getNodes() const
{
  return nodes_;
}
// =============================================================================
Teuchos::ArrayRCP<Point>
VIO::Mesh::Mesh::
getNodesNonConst()
{
  return nodes_;
}
// =============================================================================
void
VIO::Mesh::Mesh::
setNodes( const Teuchos::ArrayRCP<Point> nodes )
{
  nodes_ = nodes;
  return;
}
// =============================================================================
void
VIO::Mesh::Mesh::
setBoundaryNodes( const Teuchos::ArrayRCP<bool> isBoundaryNode )
{
  isBoundaryNode_ = isBoundaryNode;
  return; 
}
// =============================================================================
const Teuchos::ArrayRCP<const bool>
VIO::Mesh::Mesh::
getBoundaryNodes() const
{
  return isBoundaryNode_;
}
// =============================================================================
const Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> >
VIO::Mesh::Mesh::
getElems() const
{
  return elems_;
}
// =============================================================================
Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> >
VIO::Mesh::Mesh::
getElemsNonConst()
{
  return elems_;
}
// =============================================================================
void
VIO::Mesh::Mesh::
setElems( const Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems )
{
  elems_ = elems;
  return;
}
// =============================================================================
void
VIO::Mesh::Mesh::
setElemTypes( const Teuchos::ArrayRCP<Mesh::ElementType> elemTypes )
{
  elemTypes_ = elemTypes;
  return;
}
// =============================================================================
const Teuchos::ArrayRCP<const VIO::Mesh::Mesh::ElementType>
VIO::Mesh::Mesh::
getElemTypes() const
{
  return elemTypes_;
}
// =============================================================================
Teuchos::RCP<DoubleVector>
VIO::Mesh::Mesh::
getControlVolumes() const
{
  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  return controlVolumes_;
}
// =============================================================================
Teuchos::ArrayRCP<Point>
VIO::Mesh::Mesh::
getCircumcenters() const
{
  TEUCHOS_ASSERT( !circumcenters_.is_null() );
  return circumcenters_;
}
// =============================================================================
Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
VIO::Mesh::Mesh::
getEdgeLengths() const
{
  return edgeLengths_;
}
// =============================================================================
Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
VIO::Mesh::Mesh::
getCoedgeLengths() const
{
  return coedgeLengths_;
}
// =============================================================================
double
VIO::Mesh::Mesh::
getDomainArea() const
{
  return area_; 
}
// =============================================================================
void
VIO::Mesh::Mesh::
computeDomainArea_()
{
  TEUCHOS_ASSERT( !controlVolumes_.is_null() );

  // sum over the control volumes
  area_ = controlVolumes_->norm1();

  return;
}
// =============================================================================
void
VIO::Mesh::Mesh::
computeFvmEntities()
{ 
  // Compute the volume of the (Voronoi) control cells for each point.
  TEUCHOS_ASSERT( !elems_.is_null() );
  int numElems = elems_.size();
  
  TEUCHOS_ASSERT( !nodes_.is_null() );
  int numNodes = nodes_.size();
  Teuchos::RCP<Tpetra::Map<ORD> > map =
      Teuchos::rcp( new Tpetra::Map<ORD>( numNodes, 0, comm_ ) );

  controlVolumes_ = Teuchos::rcp( new DoubleVector( map ) );
  
  Teuchos::ArrayRCP<double> controlVolumesView = controlVolumes_->get1dViewNonConst();
  
  circumcenters_ = Teuchos::ArrayRCP<Point>( numElems );
  edgeLengths_   = Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( numElems );
  coedgeLengths_ = Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( numElems );
  
  // Run over the elements and calculate their contributions to the
  // control volumes.
  for ( int k=0; k<numElems; k++ )
  {
    Teuchos::ArrayRCP<ORD> & elem = elems_[k];
    
    TEST_FOR_EXCEPTION( elems_[k].size() != 3,
                        std::runtime_error,
                        "Control volumes can only be constructed consistently with triangular elements." );
    
    // extract the points for short-hand notation
    Point & x0 = nodes_[ elem[0] ];
    Point & x1 = nodes_[ elem[1] ];
    Point & x2 = nodes_[ elem[2] ];
    
    // compute the circumcenter    
    circumcenters_[k] = this->computeCircumcenter_( x0, x1, x2 );    
    Point & cc = circumcenters_[k];
    
    // midpoints of the adjacent edges
    Point mp0 = this->add_( 0.5, x0, 0.5, x1 );
    Point mp1 = this->add_( 0.5, x1, 0.5, x2 );
    Point mp2 = this->add_( 0.5, x2, 0.5, x0 );
    
    double localContributions[3];
    localContributions[0] = this->getTriangleArea_( x0, cc, mp0 )
                          + this->getTriangleArea_( x0, cc, mp2 );
    localContributions[1] = this->getTriangleArea_( x1, cc, mp1 )
                          + this->getTriangleArea_( x1, cc, mp0 );
    localContributions[2] = this->getTriangleArea_( x2, cc, mp2 )
                          + this->getTriangleArea_( x2, cc, mp1 );
                          
    controlVolumesView[ elem[0] ] += localContributions[0];
    controlVolumesView[ elem[1] ] += localContributions[1];
    controlVolumesView[ elem[2] ] += localContributions[2];
                                 
    // add the lengths of the edges
    edgeLengths_[k] = Teuchos::ArrayRCP<double>( 3 );
    edgeLengths_[k][0] = this->norm2_( this->add_( 1.0, x1, -1.0, x0 ) );
    edgeLengths_[k][1] = this->norm2_( this->add_( 1.0, x2, -1.0, x1 ) );
    edgeLengths_[k][2] = this->norm2_( this->add_( 1.0, x0, -1.0, x2 ) );
    
    // add the length of the part of the coedge in the current element
    coedgeLengths_[k] = Teuchos::ArrayRCP<double>( 3 );
    coedgeLengths_[k][0] = this->norm2_( this->add_( 1.0, mp0, -1.0, cc ) );
    coedgeLengths_[k][1] = this->norm2_( this->add_( 1.0, mp1, -1.0, cc ) );
    coedgeLengths_[k][2] = this->norm2_( this->add_( 1.0, mp2, -1.0, cc ) );

  }
  
  // TODO move this to another spot
  this->computeDomainArea_();

  return;
}
// =============================================================================
Point
VIO::Mesh::Mesh::
add_( double alpha, const Point & x,
      double beta,  const Point & y
    ) const
{
  Point z;
  for ( int k=0; k<z.size(); k++ )
      z[k] = alpha*x[k] + beta*y[k];
  
  return z;
}
// =============================================================================
double
VIO::Mesh::Mesh::
getTriangleArea_( const Point & x0,
                  const Point & x1,
                  const Point & x2
                ) const
{
    return 0.5 * this->norm2_( this->cross_( this->add_( 1.0, x1, -1.0, x0),
                                             this->add_( 1.0, x2, -1.0, x0 )
                                           )
                             );
}
// =============================================================================
Point
VIO::Mesh::Mesh::
computeCircumcenter_( const Point & x0, const Point & x1, const Point & x2
                    ) const
{ 
  Point cc;
  
  double omega = 2.0 * pow( this->norm2_( this->cross_( this->add_( 1.0, x0, -1.0, x1 ),
                                                        this->add_( 1.0, x1, -1.0, x2 ) )
                                        ), 2 );
  
  // don't divide by 0
  TEUCHOS_ASSERT( fabs(omega) > 1.0e-10 );
  
  double alpha = this->dot_( this->add_( 1.0, x1, -1.0, x2 ), this->add_( 1.0, x1, -1.0, x2 ) )
               * this->dot_( this->add_( 1.0, x0, -1.0, x1 ), this->add_( 1.0, x0, -1.0, x2 ) )
               / omega;
  double beta  = this->dot_( this->add_( 1.0, x2, -1.0, x0 ), this->add_( 1.0, x2, -1.0, x0 ) )
               * this->dot_( this->add_( 1.0, x1, -1.0, x2 ), this->add_( 1.0, x1, -1.0, x0 ) )
               / omega;
  double gamma = this->dot_( this->add_( 1.0, x0, -1.0, x1 ), this->add_( 1.0, x0, -1.0, x1 ) )
               * this->dot_( this->add_( 1.0, x2, -1.0, x0 ), this->add_( 1.0, x2, -1.0, x1 ) )
               / omega;
               
  cc = this->add_( alpha, x0, beta, x1 );
  cc = this->add_( 1.0, cc, gamma, x2 );
  
  return cc;
}
// =============================================================================
double
VIO::Mesh::Mesh::
dot_( const Point & v, const Point & w
    ) const
{
   double sum = 0.0;
   for ( int k=0; k<v.size(); k++ )
       sum += v[k] * w[k];
   return sum;
}
// =============================================================================
Point
VIO::Mesh::Mesh::
cross_( const Point & v, const Point & w
      ) const
{ 
  Point z;
  
  z[0] = v[1]*w[2] - v[2]*w[1];
  z[1] = v[2]*w[0] - v[0]*w[2];
  z[2] = v[0]*w[1] - v[1]*w[0];
  
  return z;
}
// =============================================================================
double
VIO::Mesh::Mesh::
norm2_( const Point & x
      ) const
{
  double sum = 0.0;
  for ( int k=0; k<x.size(); k++ )
      sum += x[k]*x[k];
  return sqrt( sum );
}
// =============================================================================