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

#include "VIO_TpetraMesh_Mesh.h"
#include <Epetra_Map.h>

// =============================================================================
VIO::TpetraMesh::Mesh::
Mesh( const Teuchos::RCP<const Teuchos::Comm<int> > & comm ):
    comm_( comm ),
    elems_( Teuchos::null ),
    elemTypes_( Teuchos::null ),
    nodes_( Teuchos::null ),
    isBoundaryNode_( Teuchos::null ),
    controlVolumes_( Teuchos::null ),
    area_( 0.0 ),
    scaling_( Teuchos::tuple( 1.0, 1.0, 1.0 ) ),
    fvmEntitiesUpToDate_( false )
{
}
// =============================================================================
VIO::TpetraMesh::Mesh::
~Mesh()
{
}
// =============================================================================
unsigned int
VIO::TpetraMesh::Mesh::
getNumNodes() const
{
  return nodes_.size();
}
// =============================================================================
const Teuchos::ArrayRCP<Point>
VIO::TpetraMesh::Mesh::
getNodes() const
{
  return nodes_;
}
// =============================================================================
Teuchos::ArrayRCP<Point>
VIO::TpetraMesh::Mesh::
getNodesNonConst()
{
  return nodes_;
}
// =============================================================================
void
VIO::TpetraMesh::Mesh::
setNodes( const Teuchos::ArrayRCP<Point> nodes )
{
  nodes_ = nodes;
  return;
}
// =============================================================================
void
VIO::TpetraMesh::Mesh::
setBoundaryNodes( const Teuchos::ArrayRCP<bool> isBoundaryNode )
{
  isBoundaryNode_ = isBoundaryNode;
  return;
}
// =============================================================================
const Teuchos::ArrayRCP<const bool>
VIO::TpetraMesh::Mesh::
getBoundaryNodes() const
{
  return isBoundaryNode_;
}
// =============================================================================
const Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> >
VIO::TpetraMesh::Mesh::
getElems() const
{
  return elems_;
}
// =============================================================================
Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> >
VIO::TpetraMesh::Mesh::
getElemsNonConst()
{
  return elems_;
}
// =============================================================================
void
VIO::TpetraMesh::Mesh::
setElems( const Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems )
{
  elems_ = elems;
  return;
}
// =============================================================================
void
VIO::TpetraMesh::Mesh::
setElemTypes( const Teuchos::ArrayRCP<Mesh::ElementType> elemTypes )
{
  elemTypes_ = elemTypes;
  return;
}
// =============================================================================
const Teuchos::ArrayRCP<const VIO::TpetraMesh::Mesh::ElementType>
VIO::TpetraMesh::Mesh::
getElemTypes() const
{
  return elemTypes_;
}
// =============================================================================
Teuchos::RCP<DoubleVector>
VIO::TpetraMesh::Mesh::
getControlVolumes() const
{
  if ( !fvmEntitiesUpToDate_ )
      this->computeFvmEntities_();

  TEUCHOS_ASSERT( !controlVolumes_.is_null() );

  return controlVolumes_;
}
// =============================================================================
Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
VIO::TpetraMesh::Mesh::
getEdgeLengths() const
{
  if ( !fvmEntitiesUpToDate_ )
      this->computeFvmEntities_();

  return edgeLengths_;
}
// =============================================================================
Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
VIO::TpetraMesh::Mesh::
getCoedgeLengths() const
{
  if ( !fvmEntitiesUpToDate_ )
      this->computeFvmEntities_();

  return coedgeLengths_;
}
// =============================================================================
double
VIO::TpetraMesh::Mesh::
getDomainArea() const
{
  if ( !fvmEntitiesUpToDate_ )
      this->computeFvmEntities_();

  return area_;
}
// =============================================================================
void
VIO::TpetraMesh::Mesh::
scale( double alpha )
{
   this->scale( Teuchos::tuple( alpha, alpha, alpha ) );
   return;
}
// =============================================================================
void
VIO::TpetraMesh::Mesh::
scale( const Teuchos::Tuple<double,3> & newScaling )
{
   // adapt the position of the nodes
   for ( int i=0; i<3; i++ )
   {
       if ( newScaling[i] != scaling_[i] )
       {
          double ratio = newScaling[i] / scaling_[i];
          for ( int k=0; k<nodes_.size(); k++ )
              nodes_[k][i] *= ratio;

          // store the new scaling
          scaling_[i] = newScaling[i];

          // make sure the FVM entities get updated properly
          fvmEntitiesUpToDate_ = false;
       }
   }

   return;
}
// =============================================================================
double
VIO::TpetraMesh::Mesh::
computeDomainArea_() const
{
  // break it down into a non-overlapping map
  Teuchos::RCP<Tpetra::Map<ORD> > nonoverlapMap =
      Teuchos::rcp( new Tpetra::Map<ORD>( controlVolumes_->getGlobalLength(), 0, comm_ ) );
  Teuchos::RCP<DoubleVector> tmp =
      Teuchos::rcp( new DoubleVector( nonoverlapMap ) );
  // merge the stuff into a non-overlapping vector
  Tpetra::Export<ORD> exporter( controlVolumes_->getMap(),
                                nonoverlapMap
                              );
  tmp->doExport( *controlVolumes_, exporter, Tpetra::REPLACE );

  // sum over all entries
  return tmp->norm1();
}
// =============================================================================
void
VIO::TpetraMesh::Mesh::
computeFvmEntities_() const
{
  // Compute the volume of the (Voronoi) control cells for each point.
  TEUCHOS_ASSERT( !elems_.is_null() );
  int numElems = elems_.size();

  Teuchos::RCP<Tpetra::Map<ORD> > map = this->getElemsToNodesMap_();
  controlVolumes_ = Teuchos::rcp( new DoubleVector( map ) );
  edgeLengths_    = Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( numElems );
  coedgeLengths_  = Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( numElems );

//   Teuchos::ArrayRCP<double> controlVolumesView = controlVolumes_->get1dViewNonConst();

  // Run over the elements and calculate their contributions to the
  // control volumes.
  Teuchos::ArrayRCP<double> cvView = controlVolumes_->get1dViewNonConst();
  for ( int k=0; k<numElems; k++ )
  {
    Teuchos::ArrayRCP<ORD> & elem = elems_[k];

    TEST_FOR_EXCEPTION( elems_[k].size() != 3,
                        std::runtime_error,
                        "Control volumes can only be constructed consistently with triangular elements."
                      );

    // compute the circumcenter
    Point cc = this->computeCircumcenter_( nodes_[elem[0]],
                                           nodes_[elem[1]],
                                           nodes_[elem[2]]
                                         );

    edgeLengths_[k]   = Teuchos::ArrayRCP<double>( 3 );
    coedgeLengths_[k] = Teuchos::ArrayRCP<double>( 3 );
    // iterate over the edges
    for ( int l=0; l<3; l++ )
    {
        int i0 = elem[ l ];
        int i1 = elem[ (l+1)%3 ];

        Point & x0 = nodes_[i0];
        Point & x1 = nodes_[i1];

        // edge midpoint
        Point mp = this->add_( 0.5, x0, 0.5, x1 );

        cvView[ map->getLocalElement(i0) ] += this->getTriangleArea_( x0, cc, mp );
        cvView[ map->getLocalElement(i1) ] += this->getTriangleArea_( x1, cc, mp );

        coedgeLengths_[k][l] = this->norm2_( this->add_( 1.0, mp, -1.0, cc ) );
        edgeLengths_[k][l]   = this->norm2_( this->add_( 1.0, x1, -1.0, x0 ) );
    }
  }

  // sum up the overlapping entries and make sure they're
  // available on all processes
  this->sumInOverlapMap_( controlVolumes_ );

  // TODO move this to another spot
  area_ = this->computeDomainArea_();

  fvmEntitiesUpToDate_ = true;

  return;
}
// =============================================================================
void
VIO::TpetraMesh::Mesh::
sumInOverlapMap_( Teuchos::RCP<DoubleVector> x ) const
{
  Teuchos::RCP<Tpetra::Map<ORD> > nonoverlapMap =
      Teuchos::rcp( new Tpetra::Map<ORD>( x->getGlobalLength(), 0, comm_ ) );
  Teuchos::RCP<DoubleVector> tmp =
      Teuchos::rcp( new DoubleVector( nonoverlapMap ) );

  // merge the stuff into a non-overlapping vector
  Tpetra::Export<ORD> exporter( x->getMap(),
                                nonoverlapMap
                              );
  tmp->doExport( *x, exporter, Tpetra::ADD );

  // map it back out to x
  x->doImport( *tmp, exporter, Tpetra::REPLACE );

  return;
}
// =============================================================================
Teuchos::RCP<Tpetra::Map<ORD> >
VIO::TpetraMesh::Mesh::
getElemsToNodesMap_() const
{
  // create list of elements that need to be accessible from this process


  // Make sure that *all entries that belong to any of the elements in
  // this core are accessible.
  // First mark all the nodes that need to be accessible:
  TEUCHOS_ASSERT( !elems_.is_null() );
  TEUCHOS_ASSERT( !nodes_.is_null() );
  int numNodes = nodes_.size();
  Teuchos::Array<bool> mustBeAccessible( numNodes );
  for ( int k=0; k<elems_.size(); k++ )
      for ( int l=0; l<elems_[k].size(); l++ )
          mustBeAccessible[ elems_[k][l] ] = true;
  // now create the list
  Teuchos::Array<ORD> entryList;
  for ( int k=0; k<numNodes; k++ )
      if ( mustBeAccessible[k] )
          entryList.append( k );

  Teuchos::RCP<Tpetra::Map<ORD> > map =
      Teuchos::rcp( new Tpetra::Map<ORD>( Teuchos::OrdinalTraits<ORD>::invalid(), entryList(), 0, comm_ )
                  );

  return map;
}
// =============================================================================
Point
VIO::TpetraMesh::Mesh::
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
VIO::TpetraMesh::Mesh::
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
VIO::TpetraMesh::Mesh::
computeCircumcenter_( const Point & x0, const Point & x1, const Point & x2
                    ) const
{
  Point cc;

  double omega = 2.0 * pow( this->norm2_( this->cross_( this->add_( 1.0, x0, -1.0, x1 ),
                                                        this->add_( 1.0, x1, -1.0, x2 ) )
                                        ), 2 );

  // don't divide by 0
  TEST_FOR_EXCEPTION( fabs(omega) < 1.0e-10,
                      std::runtime_error,
                      "It seems that the vectors \n\n"
                      << "   " << x0 << "\n"
                      << "   " << x1 << "\n"
                      << "   " << x2 << "\n"
                      << "\ndo not form a proper triangle. Abort."
                      << std::endl
                    );

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
VIO::TpetraMesh::Mesh::
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
VIO::TpetraMesh::Mesh::
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
VIO::TpetraMesh::Mesh::
norm2_( const Point & x
      ) const
{
  double sum = 0.0;
  for ( int k=0; k<x.size(); k++ )
      sum += x[k]*x[k];
  return sqrt( sum );
}
// =============================================================================
