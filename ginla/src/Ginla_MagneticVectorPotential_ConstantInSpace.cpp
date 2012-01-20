// @HEADER
//
//    Query routines for the magnetic vector potential.
//    Copyright (C) 2012  Nico Schl\"omer
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

#include "Ginla_MagneticVectorPotential_ConstantInSpace.hpp"
#include "Ginla_StkMesh.hpp"

#include <Epetra_Vector.h>

namespace Ginla {
namespace MagneticVectorPotential {
// ============================================================================
ConstantInSpace::
ConstantInSpace( const Teuchos::RCP<Ginla::StkMesh> &mesh,
                 const Teuchos::RCP<DoubleVector> &b,
                 double mu,
                 double theta,
                 const Teuchos::RCP<DoubleVector> &u
                 ) :
  mesh_( mesh ),
  b_( b ),
  rotatedB_( *b ),
  dRotatedBDTheta_( DoubleVector( 3 ) ),
  mu_( mu ),
  theta_( theta ),
  u_( u ),
  edgeCache_( Teuchos::ArrayRCP<DoubleVector>() ),
  edgeCacheUptodate_( false ),
  edgeCacheFallback_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> >() ),
  edgeCacheFallbackUptodate_( false )
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !b_.is_null() );
#endif
  TEST_FOR_EXCEPT_MSG( b_->dot( *b_ ) != 1.0,
                       "Magnetic field vector not normalized: "
                       << "<b,b> = " << b->dot( *b ) << "."
                       << std::endl
                       );
  if ( !u_.is_null() )
  {
    TEST_FOR_EXCEPT_MSG( u_->dot( *u_ ) != 1.0,
                         "Rotation vector not normalized: "
                         << "<u,u> = " << u_->dot( *u_ ) << "."
                         << std::endl
                         );
    double sinTheta, cosTheta;
    sincos( theta_, &sinTheta, &cosTheta );
    rotatedB_ = this->rotate_( *b_, *u_, sinTheta, cosTheta );
    dRotatedBDTheta_ = this->dRotateDTheta_( *b_, *u_, sinTheta, cosTheta );
  }


  if ( mesh->supportsEdges() )
  {
    edgeCache_ =
      Teuchos::ArrayRCP<DoubleVector>( mesh_->getOverlapEdges().size() );
  }
  else
  {
    edgeCacheFallback_ =
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> >( mesh_->getOwnedCells(
                                                             ).size() );
  }
  return;
}
// ============================================================================
ConstantInSpace::
~ConstantInSpace()
{
}
// ============================================================================
DoubleVector
ConstantInSpace::
rotate_( const DoubleVector &v,
         const DoubleVector &u,
         const double sinTheta,
         const double cosTheta
         ) const
{
  // Rotate a vector \c v by the angle \c theta in the plane perpendicular
  // to the axis given by \c u.
  // Refer to
  // http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

  DoubleVector r = v;
  if ( sinTheta != 0.0 )
  {
    // cos(theta) * I * v
    r *= cosTheta;

    // + sin(theta) u\cross v
    DoubleVector tmp = this->crossProduct_( u, v );
    tmp *= sinTheta;
    r += tmp;

    // + (1-cos(theta)) (u*u^T) * v
    tmp = u;
    tmp *= (1.0-cosTheta) * u.dot( v );

    r += tmp;
  }

  return r;
}
// ============================================================================
DoubleVector
ConstantInSpace::
dRotateDTheta_( const DoubleVector &v,
                const DoubleVector &u,
                const double sinTheta,
                const double cosTheta
                ) const
{
  DoubleVector r = v;

  // cos(theta) * I * v
  r *= -sinTheta;
  // + sin(theta) u\cross v
  DoubleVector tmp = this->crossProduct_( u, v );
  tmp *= cosTheta;
  r += tmp;
  // + (1-cos(theta)) (u*u^T) * v
  tmp = u;
  tmp *= sinTheta * u.dot( v );
  r += tmp;

  return r;
}
// ============================================================================
DoubleVector
ConstantInSpace::
crossProduct_( const DoubleVector u,
               const DoubleVector v
               ) const
{
  DoubleVector uXv( 3 );
  uXv[0] = u[1]*v[2] - u[2]*v[1];
  uXv[1] = u[2]*v[0] - u[0]*v[2];
  uXv[2] = u[0]*v[1] - u[1]*v[0];
  return uXv;
}
// ============================================================================
void
ConstantInSpace::
setParameters( const LOCA::ParameterVector &p )
{
  if (p.isParameter( "mu" ))
    mu_ = p.getValue( "mu" );

  if (p.isParameter( "theta" ))
  {
    theta_ = p.getValue( "theta" );
    double sinTheta, cosTheta;
    sincos( theta_, &sinTheta, &cosTheta );
    rotatedB_ = this->rotate_( *b_, *u_, sinTheta, cosTheta );
    dRotatedBDTheta_ = this->dRotateDTheta_( *b_, *u_, sinTheta, cosTheta );
  }

  return;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
ConstantInSpace::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
    Teuchos::rcp( new LOCA::ParameterVector() );

  p->addParameter( "mu", mu_ );
  p->addParameter( "theta", theta_ );

  return p;
}
// ============================================================================
double
ConstantInSpace::
getAEdgeMidpointProjection( const unsigned int edgeIndex
                            ) const
{
  // A vector potential associated with the constant magnetic field RB is
  //
  //     A(X) = 0.5 * (RB x X).
  //
  // Projecting it onto an edge e yields
  //     e * A(X) = e * (0.5*(RB x X))
  //              = RB * (0.5*(X x e)),
  //
  // taking this at the edge midpoint Em of the edge e gives
  //     e * A(Em) = RB * (0.5*(Em x e)).
  //
  // The vector 0.5*(Em x e) can be cached for the mesh.
  // This saves caching e and the edge midpoint separately and also avoids
  // computing the cross-products more than once if B changes.
  if ( !edgeCacheUptodate_ )
    this->initializeEdgeCache_();

  return mu_ * rotatedB_.dot( edgeCache_[edgeIndex] );
}
// ============================================================================
double
ConstantInSpace::
getdAdMuEdgeMidpointProjection( const unsigned int edgeIndex
                                ) const
{
  if ( !edgeCacheUptodate_ )
    this->initializeEdgeCache_();

  return rotatedB_.dot( edgeCache_[edgeIndex] );
}
// ============================================================================
double
ConstantInSpace::
getdAdThetaEdgeMidpointProjection( const unsigned int edgeIndex
                                   ) const
{
  if ( !edgeCacheUptodate_ )
    this->initializeEdgeCache_();

  return mu_ * dRotatedBDTheta_.dot( edgeCache_[edgeIndex] );
}
// ============================================================================
void
ConstantInSpace::
initializeEdgeCache_() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  std::vector<stk::mesh::Entity*> edges = mesh_->getOverlapEdges();

  // Loop over all edges and create the cache.
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
    // Get the two end points.
    const stk::mesh::PairIterRelation endPoints =
      edges[k]->relations( mesh_->getMetaData()->node_rank() );

    // extract the nodal coordinates
    Teuchos::ArrayRCP<DoubleVector> localNodeCoords =
      mesh_->getNodeCoordinates( endPoints );
#ifdef _DEBUG_
    TEUCHOS_ASSERT_EQUALITY( localNodeCoords.size(), 2 );
#endif
    DoubleVector edge = localNodeCoords[1];
    edge -= localNodeCoords[0];
    DoubleVector edgeMidpoint = localNodeCoords[1];
    edgeMidpoint += localNodeCoords[0];
    TEUCHOS_ASSERT_EQUALITY( 0, edgeMidpoint.scale( 0.5 ));

    // Compute the cross product.
    edgeCache_[k] = this->crossProduct_( edgeMidpoint, edge );
    edgeCache_[k].scale( 0.5 );
  }

  edgeCacheUptodate_ = true;

  return;
}
// ============================================================================
double
ConstantInSpace::
getAEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                    const unsigned int edgeIndex
                                    ) const
{
  if ( !edgeCacheFallbackUptodate_ )
    this->initializeEdgeCacheFallback_();

  return mu_ * rotatedB_.dot( edgeCacheFallback_[cellIndex][edgeIndex] );
}
// ============================================================================
double
ConstantInSpace::
getdAdMuEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                        const unsigned int edgeIndex
                                        ) const
{
  if ( !edgeCacheFallbackUptodate_ )
    this->initializeEdgeCacheFallback_();

  return rotatedB_.dot( edgeCacheFallback_[cellIndex][edgeIndex] );
}
// ============================================================================
double
ConstantInSpace::
getdAdThetaEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                           const unsigned int edgeIndex
                                           ) const
{
  if ( !edgeCacheFallbackUptodate_ )
    this->initializeEdgeCacheFallback_();

  return mu_ * dRotatedBDTheta_.dot( edgeCacheFallback_[cellIndex][edgeIndex] );
}
// ============================================================================
void
ConstantInSpace::
initializeEdgeCacheFallback_() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !edgeCacheFallback_.is_null() );
#endif

  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // Loop over all edges and create the cache.
  // To this end, loop over all cells and the edges within the cell.
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
    // get the nodes local to the cell
    stk::mesh::PairIterRelation localNodes = (*cells[k]).relations(
        mesh_->getMetaData()->node_rank() );

    unsigned int numLocalNodes = localNodes.size();
    unsigned int cellDimension = mesh_->getCellDimension( numLocalNodes );
    // extract the nodal coordinates
    Teuchos::ArrayRCP<DoubleVector> localNodeCoords =
      mesh_->getNodeCoordinates( localNodes );

    edgeCacheFallback_[k] =
      Teuchos::ArrayRCP<DoubleVector>( mesh_->getNumEdgesPerCell( cellDimension ) );

    // In a simplex, the edges are exactly the connection between each pair
    // of nodes. Hence, loop over pairs of nodes.
    unsigned int edgeIndex = 0;
    for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
    {
      for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
      {
        DoubleVector edge = localNodeCoords[1];
        edge -= localNodeCoords[0];
        DoubleVector edgeMidpoint = localNodeCoords[1];
        edgeMidpoint += localNodeCoords[0];
        TEUCHOS_ASSERT_EQUALITY( 0, edgeMidpoint.scale( 0.5 ));

        // Compute the cross product.
        edgeCacheFallback_[k][edgeIndex] =
          this->crossProduct_( edgeMidpoint, edge );
        edgeCacheFallback_[k][edgeIndex].scale( 0.5 );

        edgeIndex++;
      }
    }
  }

  edgeCacheFallbackUptodate_ = true;
  return;
}
// ============================================================================
} // namespace MagneticVectorPotential
} // namespace Ginla
