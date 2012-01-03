// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2011  Nico Schl\"omer
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

#include "Ginla_MagneticVectorPotential.hpp"
#include "Ginla_StkMesh.hpp"

#include <Epetra_Vector.h>

namespace Ginla {
// ============================================================================
MagneticVectorPotential::
MagneticVectorPotential( const Teuchos::RCP<Ginla::StkMesh>           & mesh,
                         const Teuchos::RCP<const Epetra_MultiVector> & mvp,
                         double mu,
                         double theta,
                         const DoubleVector u
                       ):
  mesh_( mesh ),
  mvp_( mvp ),
  mu_( mu ),
  theta_( theta ),
  sinTheta_( 0.0 ),
  cosTheta_( 1.0 ),
  u_( u ),
  mvpEdgeMidpoint_( Teuchos::ArrayRCP<DoubleVector>() ),
  edges_( Teuchos::ArrayRCP<DoubleVector>() ),
  mvpEdgeMidpointUpToDate_( false ),
  mvpEdgeMidpointFallback_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> >() ),
  edgesFallback_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> >() ),
  mvpEdgeMidpointFallbackUpToDate_( false )
{
    if ( theta_ != 0.0 )
        sincos( theta_, &sinTheta_, &cosTheta_ );

    // Normalize u.
    double uDotu = u_.dot(u_);
    if ( uDotu != 0.0 && uDotu != 1.0 )
        u_ *= 1.0 / sqrt(uDotu);

#ifdef _DEBUG_
    TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
    try
    {
        mvpEdgeMidpoint_ = Teuchos::ArrayRCP<DoubleVector>( mesh_->getOverlapEdges().size() );
        edges_ = Teuchos::ArrayRCP<DoubleVector>( mesh_->getOverlapEdges().size() );
    }
    catch( ... )
    {
        mvpEdgeMidpointFallback_ =
            Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> >( mesh_->getOwnedCells().size() );
        edgesFallback_ =
            Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> >( mesh_->getOwnedCells().size() );
    }
    return;
}
// ============================================================================
MagneticVectorPotential::
~MagneticVectorPotential()
{
}
// ============================================================================
void
MagneticVectorPotential::
setParameters( const LOCA::ParameterVector & p )
{
    if (p.isParameter( "mu" ))
        mu_ = p.getValue ( "mu" );

    if (p.isParameter( "theta" ))
    {
        theta_ = p.getValue ( "theta" );
        sincos( theta_, &sinTheta_, &cosTheta_ );
    }

    return;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
MagneticVectorPotential::
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
MagneticVectorPotential::
getAEdgeMidpointProjection( const unsigned int edgeIndex
                          ) const
{
  if ( !mvpEdgeMidpointUpToDate_ )
      this->initializeMvpEdgeMidpointCache_();

  // perform rotation
  DoubleVector x = this->rotate_( mvpEdgeMidpoint_[edgeIndex],
                                  u_, sinTheta_, cosTheta_
                                );

  return mu_ * edges_[edgeIndex].dot( x );
}
// ============================================================================
DoubleVector
MagneticVectorPotential::
rotate_( const DoubleVector & v,
         const DoubleVector & u,
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
      DoubleVector tmp = this->crossProduct_(u,v);
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
MagneticVectorPotential::
crossProduct_( const DoubleVector u,
               const DoubleVector v
             ) const
{
    DoubleVector uXv;
    uXv[0] = u[1]*v[2] - u[2]*v[1];
    uXv[1] = u[2]*v[0] - u[0]*v[2];
    uXv[2] = u[0]*v[1] - u[1]*v[0];
    return uXv;
}
// ============================================================================
double
MagneticVectorPotential::
getdAdMuEdgeMidpointProjection( const unsigned int edgeIndex
                              ) const
{
  if ( !mvpEdgeMidpointUpToDate_ )
      this->initializeMvpEdgeMidpointCache_();

  // perform rotation
  DoubleVector x = this->rotate_( mvpEdgeMidpoint_[edgeIndex],
                                  u_, sinTheta_, cosTheta_
                                );

  return edges_[edgeIndex].dot( x );
}
// ============================================================================
double
MagneticVectorPotential::
getdAdThetaEdgeMidpointProjection( const unsigned int edgeIndex
                                  ) const
{
  if ( !mvpEdgeMidpointUpToDate_ )
      this->initializeMvpEdgeMidpointCache_();

  // perform rotation
  DoubleVector x = this->rotate_( mvpEdgeMidpoint_[edgeIndex],
                                  u_, cosTheta_, -sinTheta_
                                );

  return mu_ * edges_[edgeIndex].dot( x );
}
// ============================================================================
void
MagneticVectorPotential::
initializeMvpEdgeMidpointCache_() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !mvp_.is_null() );
#endif
  std::vector<stk::mesh::Entity*> edges = mesh_->getOverlapEdges();

  // Loop over all edges and create the cache.
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
      // Get the two end points.
      const stk::mesh::PairIterRelation endPoints =
          edges[k]->relations( mesh_->getMetaData()->node_rank() );

      // get the local ids
      Teuchos::Tuple<int,2> lid;
      lid[0] = mvp_->Map().LID( (*endPoints[0].entity()).identifier() - 1 );
      lid[1] = mvp_->Map().LID( (*endPoints[1].entity()).identifier() - 1 );
#ifdef _DEBUG_
      TEST_FOR_EXCEPT_MSG( lid[0] < 0,
                          "The global index " << (*endPoints[0].entity()).identifier() - 1
                          << " does not seem to be present on this node." );
      TEST_FOR_EXCEPT_MSG( lid[1] < 0,
                          "The global index " << (*endPoints[1].entity()).identifier() - 1
                          << " does not seem to be present on this node." );
#endif

      // Approximate the value at the midpoint of the edge
      // by the average of the values at the adjacent nodes.
      mvpEdgeMidpoint_[k] = DoubleVector(3);
      for (int i=0; i<3; i++ )
          mvpEdgeMidpoint_[k][i] = 0.5 * ( (*(*mvp_)(i))[lid[0]]
                                         + (*(*mvp_)(i))[lid[1]] );

      // extract the nodal coordinates
      Teuchos::ArrayRCP<DoubleVector> localNodes =
          mesh_->getNodeCoordinates( endPoints );
#ifdef _DEBUG_
      TEUCHOS_ASSERT_EQUALITY( localNodes.size(), 2 );
#endif
      edges_[k] = localNodes[1];
      edges_[k] -= localNodes[0];
  }

  mvpEdgeMidpointUpToDate_ = true;
  return;
}
// ============================================================================
double
MagneticVectorPotential::
getAEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                    const unsigned int edgeIndex
                                  ) const
{
  if ( !mvpEdgeMidpointFallbackUpToDate_ )
      this->initializeMvpEdgeMidpointFallback_();

  // perform rotation
  DoubleVector x = this->rotate_( mvpEdgeMidpointFallback_[cellIndex][edgeIndex],
                                  u_, sinTheta_, cosTheta_
                                );

  return mu_ * edgesFallback_[cellIndex][edgeIndex].dot( x );
}
// ============================================================================
double
MagneticVectorPotential::
getdAdMuEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                        const unsigned int edgeIndex
                                      ) const
{
  if ( !mvpEdgeMidpointFallbackUpToDate_ )
      this->initializeMvpEdgeMidpointFallback_();

  // perform rotation
  DoubleVector x = this->rotate_( mvpEdgeMidpointFallback_[cellIndex][edgeIndex],
                                  u_, sinTheta_, cosTheta_
                                );

  return edgesFallback_[cellIndex][edgeIndex].dot( x );
}
// ============================================================================
double
MagneticVectorPotential::
getdAdThetaEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                            const unsigned int edgeIndex
                                          ) const
{
  if ( !mvpEdgeMidpointFallbackUpToDate_ )
      this->initializeMvpEdgeMidpointFallback_();

  // perform rotation
  DoubleVector x = this->rotate_( mvpEdgeMidpoint_[edgeIndex],
                                  u_, cosTheta_, -sinTheta_
                                );

  return mu_ * edgesFallback_[cellIndex][edgeIndex].dot( x );
}
// ============================================================================
void
MagneticVectorPotential::
initializeMvpEdgeMidpointFallback_() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !mvp_.is_null() );
#endif
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // Loop over all edges and create the cache.
  // To this end, loop over all cells and the edges within the cell.
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation localNodes = (*cells[k]).relations( mesh_->getMetaData()->node_rank() );

      unsigned int numLocalNodes = localNodes.size();
      unsigned int cellDimension = mesh_->getCellDimension( numLocalNodes );
      // extract the nodal coordinates
      Teuchos::ArrayRCP<DoubleVector> localNodeCoords =
          mesh_->getNodeCoordinates( localNodes );

      mvpEdgeMidpointFallback_[k] =
          Teuchos::ArrayRCP<DoubleVector>( mesh_->getNumEdgesPerCell( cellDimension ) );
      edgesFallback_[k] =
          Teuchos::ArrayRCP<DoubleVector>( mesh_->getNumEdgesPerCell( cellDimension ) );

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      unsigned int edgeIndex = 0;
      Teuchos::Tuple<int,2> lid;
      for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
      {
          lid[0] = mvp_->Map().LID( (*localNodes[e0].entity()).identifier() - 1 );
#ifdef _DEBUG_
          TEST_FOR_EXCEPT_MSG( lid[0] < 0,
                               "The global index " << (*localNodes[e0].entity()).identifier() - 1
                               << " does not seem to be present on this node." );
#endif
          for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
          {
              lid[1] = mvp_->Map().LID( (*localNodes[e1].entity()).identifier() - 1 );
#ifdef _DEBUG_
              TEST_FOR_EXCEPT_MSG( lid[1] < 0,
                                   "The global index " << (*localNodes[e1].entity()).identifier() - 1
                                   << " does not seem to be present on this node." );
#endif

              // Approximate the value at the midpoint of the edge
              // by the average of the values at the adjacent nodes.
              mvpEdgeMidpointFallback_[k][edgeIndex] = DoubleVector(3);
              for (int i=0; i<3; i++ )
                  mvpEdgeMidpointFallback_[k][edgeIndex][i] = 0.5 * ( (*(*mvp_)(i))[lid[0]]
                                                                    + (*(*mvp_)(i))[lid[1]] );

              // Update edge cache.
              edgesFallback_[k][edgeIndex] = localNodeCoords[e1];
              edgesFallback_[k][edgeIndex] -= localNodeCoords[e0];

              edgeIndex++;
          }
      }
  }

  mvpEdgeMidpointFallbackUpToDate_ = true;
  return;
}
// ============================================================================
} // namespace Ginla
