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

#include "Ginla_MagneticVectorPotential_ExplicitValues.hpp"
#include "Ginla_StkMesh.hpp"

#include <Epetra_Vector.h>

namespace Ginla {
namespace MagneticVectorPotential {
// ============================================================================
ExplicitValues::
ExplicitValues( const Teuchos::RCP<Ginla::StkMesh>           & mesh,
                const Teuchos::RCP<const Epetra_MultiVector> & mvp,
                double mu
              ):
  mesh_( mesh ),
  mvp_( mvp ),
  mu_( mu ),
  mvpEdgeMidpointProjectionCache_( Teuchos::ArrayRCP<double>() ),
  mvpEdgeMidpointProjectionCacheUptodate_( false ),
  mvpEdgeMidpointProjectionCacheFallback_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >() ),
  mvpEdgeMidpointProjectionCacheFallbackUptodate_( false )
{
#ifdef _DEBUG_
    TEUCHOS_ASSERT( !mesh_.is_null() );
#endif

    if ( mesh->supportsEdges() )
    {
        mvpEdgeMidpointProjectionCache_ =
            Teuchos::ArrayRCP<double>( mesh_->getOverlapEdges().size() );
    }
    else
    {
        mvpEdgeMidpointProjectionCacheFallback_ =
            Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( mesh_->getOwnedCells().size() );
    }

    return;
}
// ============================================================================
ExplicitValues::
~ExplicitValues()
{
}
// ============================================================================
void
ExplicitValues::
setParameters( const LOCA::ParameterVector & p )
{
    if (p.isParameter( "mu" ))
        mu_ = p.getValue ( "mu" );

    return;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
ExplicitValues::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );

  p->addParameter( "mu", mu_ );

  return p;
}
// ============================================================================
double
ExplicitValues::
getAEdgeMidpointProjection( const unsigned int edgeIndex
                          ) const
{
  if ( !mvpEdgeMidpointProjectionCacheUptodate_ )
      this->initializeMvpEdgeMidpointCache_();

  return mu_ * mvpEdgeMidpointProjectionCache_[edgeIndex];
}
// ============================================================================
double
ExplicitValues::
getdAdMuEdgeMidpointProjection( const unsigned int edgeIndex
                              ) const
{
  if ( !mvpEdgeMidpointProjectionCacheUptodate_ )
      this->initializeMvpEdgeMidpointCache_();

  return mvpEdgeMidpointProjectionCache_[edgeIndex];
}
// ============================================================================
double
ExplicitValues::
getdAdThetaEdgeMidpointProjection( const unsigned int edgeIndex
                                  ) const
{
  TEST_FOR_EXCEPT_MSG( true,
                       "Parameter \"theta\" not supported. Abort." );
}
// ============================================================================
void
ExplicitValues::
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
      DoubleVector mvpEdgeMidpoint(3);
      for (int i=0; i<3; i++ )
          mvpEdgeMidpoint[i] = 0.5 * ( (*(*mvp_)(i))[lid[0]]
                                     + (*(*mvp_)(i))[lid[1]] );

      // extract the nodal coordinates
      Teuchos::ArrayRCP<DoubleVector> localNodes =
          mesh_->getNodeCoordinates( endPoints );
#ifdef _DEBUG_
      TEUCHOS_ASSERT_EQUALITY( localNodes.size(), 2 );
#endif
      DoubleVector edge = localNodes[1];
      edge -= localNodes[0];

      mvpEdgeMidpointProjectionCache_[k] = edge.dot( mvpEdgeMidpoint );
  }

  mvpEdgeMidpointProjectionCacheUptodate_ = true;
  return;
}
// ============================================================================
double
ExplicitValues::
getAEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                    const unsigned int edgeIndex
                                  ) const
{
  if ( !mvpEdgeMidpointProjectionCacheFallbackUptodate_ )
      this->initializeMvpEdgeMidpointFallback_();

  return mu_ * mvpEdgeMidpointProjectionCacheFallback_[cellIndex][edgeIndex];
}
// ============================================================================
double
ExplicitValues::
getdAdMuEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                        const unsigned int edgeIndex
                                      ) const
{
  if ( !mvpEdgeMidpointProjectionCacheFallbackUptodate_ )
      this->initializeMvpEdgeMidpointFallback_();

  return mvpEdgeMidpointProjectionCacheFallback_[cellIndex][edgeIndex];
}
// ============================================================================
double
ExplicitValues::
getdAdThetaEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                           const unsigned int edgeIndex
                                         ) const
{
  TEST_FOR_EXCEPT_MSG( true,
                       "Parameter \"theta\" not supported. Abort." );
}
// ============================================================================
void
ExplicitValues::
initializeMvpEdgeMidpointFallback_() const
{
// #ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !mvp_.is_null() );
  TEUCHOS_ASSERT( !mvpEdgeMidpointProjectionCacheFallback_.is_null() );
// #endif

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

      mvpEdgeMidpointProjectionCacheFallback_[k] =
          Teuchos::ArrayRCP<double>( mesh_->getNumEdgesPerCell( cellDimension ) );

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
              DoubleVector mvpEdgeMidpoint(3);
              for (int i=0; i<3; i++ )
                  mvpEdgeMidpoint[i] =
                      0.5 * ( (*(*mvp_)(i))[lid[0]] + (*(*mvp_)(i))[lid[1]] );

              // Update edge cache.
              DoubleVector edge = localNodeCoords[e1];
              edge -= localNodeCoords[e0];

              mvpEdgeMidpointProjectionCacheFallback_[k][edgeIndex] =
                  edge.dot( mvpEdgeMidpoint );

              edgeIndex++;
          }
      }
  }

  mvpEdgeMidpointProjectionCacheFallbackUptodate_ = true;
  return;
}
// ============================================================================
} // namespace MagneticVectorPotential
} // namespace Ginla
