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
                         double mu
                       ):
  mesh_( mesh ),
  mvp_( mvp ),
  mu_( mu ),
  edgeMidpointProjectionCache_( Teuchos::ArrayRCP<double>() ),
  edgeMidpointProjectionCacheUpToDate_( false ),
  edgeMidpointProjectionFallbackCache_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >() ),
  edgeMidpointProjectionFallbackCacheUpToDate_( false )
{
    TEUCHOS_ASSERT( !mesh_.is_null() );
    // TODO remove one
    try
    {
        edgeMidpointProjectionCache_ = Teuchos::ArrayRCP<double>( mesh_->getOverlapEdges().size() );
    }
    catch( ... )
    {
        edgeMidpointProjectionCache_ = Teuchos::ArrayRCP<double>();
    }
    edgeMidpointProjectionFallbackCache_ = Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( mesh_->getOwnedCells().size() );
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
        if ( mu_ != p.getValue ( "mu" ) )
        {
            mu_ = p.getValue ( "mu" );
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

  return p;
}
// ============================================================================
Teuchos::RCP<DoubleVector>
MagneticVectorPotential::
getA(const DoubleVector & x) const
{
  TEST_FOR_EXCEPT_MSG( true,
                               "getA(x) for general 'x' cannot be implemented as "
                               << "no general MVP info is provided to this class."
                             );
  Teuchos::RCP<DoubleVector> a;
  return a;
}
// ============================================================================
double
MagneticVectorPotential::
getAEdgeMidpointProjection( const unsigned int edgeIndex
                          ) const
{
  if ( !edgeMidpointProjectionCacheUpToDate_ )
      this->initializeEdgeMidpointProjectionCache_();

  return mu_ * edgeMidpointProjectionCache_[edgeIndex];
}
// ============================================================================
double
MagneticVectorPotential::
getdAdMuEdgeMidpointProjection( const unsigned int edgeIndex
                              ) const
{
  if ( !edgeMidpointProjectionCacheUpToDate_ )
      this->initializeEdgeMidpointProjectionCache_();

  return edgeMidpointProjectionCache_[edgeIndex];
}
// ============================================================================
void
MagneticVectorPotential::
initializeEdgeMidpointProjectionCache_() const
{
  TEUCHOS_ASSERT( !mesh_.is_null() );
  std::vector<stk::mesh::Entity*> edges = mesh_->getOverlapEdges();

  TEUCHOS_ASSERT( !mvp_.is_null() );
  // Loop over all edges and create the cache.
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
      // Get the two end points.
      stk::mesh::PairIterRelation endPoints =
          edges[k]->relations( mesh_->getMetaData()->node_rank() );

      // get the local ids
      Teuchos::Tuple<int,2> gid, lid;
      gid[0] = (*endPoints[0].entity()).identifier() - 1;
      lid[0] = mvp_->Map().LID( gid[0] );
      TEST_FOR_EXCEPT_MSG( lid[0] < 0,
                            "The global index " << gid[0]
                            << " does not seem to be present on this node." );
      gid[1] = (*endPoints[1].entity()).identifier() - 1;
      lid[1] = mvp_->Map().LID( gid[1] );
      TEST_FOR_EXCEPT_MSG( lid[1] < 0,
                            "The global index " << gid[1]
                            << " does not seem to be present on this node." );

      // Approximate the value at the midpoint of the edge
      // by the average of the values at the adjacent nodes.
      Teuchos::SerialDenseVector<int,double> a(3);
      for (int i=0; i<3; i++ )
          a[i] = 0.5 * ( (*(*mvp_)(i))[lid[0]] + (*(*mvp_)(i))[lid[1]] );

      // extract the nodal coordinates
      Teuchos::ArrayRCP<Teuchos::SerialDenseVector<int,double> > localNodes =
          mesh_->getNodeCoordinates( endPoints );
      TEUCHOS_ASSERT_EQUALITY( localNodes.size(), 2 );
      Teuchos::SerialDenseVector<int,double> edge = localNodes[1];
      edge -= localNodes[0];

      // Fill the cache.
      edgeMidpointProjectionCache_[k] = edge.dot( a );
  }

  edgeMidpointProjectionCacheUpToDate_ = true;
  return;
}
// ============================================================================
double
MagneticVectorPotential::
getAEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                    const unsigned int edgeIndex
                                  ) const
{
  if ( !edgeMidpointProjectionFallbackCacheUpToDate_ )
      this->initializeEdgeMidpointProjectionFallbackCache_();

  return mu_ * edgeMidpointProjectionFallbackCache_[cellIndex][edgeIndex];
}
// ============================================================================
double
MagneticVectorPotential::
getdAdMuEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                        const unsigned int edgeIndex
                                      ) const
{
  if ( !edgeMidpointProjectionFallbackCacheUpToDate_ )
      this->initializeEdgeMidpointProjectionFallbackCache_();

  return edgeMidpointProjectionFallbackCache_[cellIndex][edgeIndex];
}
// ============================================================================
void
MagneticVectorPotential::
initializeEdgeMidpointProjectionFallbackCache_() const
{
  TEUCHOS_ASSERT( !mesh_.is_null() );
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // Loop over all edges and create the cache.
  // To this end, loop over all cells and the edges within the cell.
  TEUCHOS_ASSERT( !mvp_.is_null() );
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation rel = (*cells[k]).relations( mesh_->getMetaData()->node_rank() );

      unsigned int numLocalNodes = rel.size();
      unsigned int cellDimension = mesh_->getCellDimension( numLocalNodes );
      // extract the nodal coordinates
      Teuchos::ArrayRCP<Teuchos::SerialDenseVector<int,double> > localNodes =
          mesh_->getNodeCoordinates( rel );

      edgeMidpointProjectionFallbackCache_[k] =
          Teuchos::ArrayRCP<double>( mesh_->getNumEdgesPerCell( cellDimension ) );

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      unsigned int edgeIndex = 0;
      Teuchos::Tuple<int,2> gid, lid;
      for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
      {
          gid[0] = (*rel[e0].entity()).identifier() - 1;
          lid[0] = mvp_->Map().LID( gid[0] );
          TEST_FOR_EXCEPT_MSG( lid[0] < 0,
                               "The global index " << gid[0]
                               << " does not seem to be present on this node." );
          for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
          {
              gid[1] = (*rel[e1].entity()).identifier() - 1;
              lid[1] = mvp_->Map().LID( gid[1] );
              TEST_FOR_EXCEPT_MSG( lid[1] < 0,
                                   "The global index " << gid[1]
                                   << " does not seem to be present on this node." );

              // Approximate the value at the midpoint of the edge
              // by the average of the values at the adjacent nodes.
              Teuchos::SerialDenseVector<int,double> a(3);
              for (int i=0; i<3; i++ )
                  a[i] = 0.5 * ( (*(*mvp_)(i))[lid[0]] + (*(*mvp_)(i))[lid[1]] );

              edgeMidpointProjectionFallbackCache_[k][edgeIndex] = 0.0;
              Teuchos::SerialDenseVector<int,double> edge = localNodes[e1];
              edge -= localNodes[e0];
              edgeMidpointProjectionFallbackCache_[k][edgeIndex++] = edge.dot( a );
          }
      }
  }

  edgeMidpointProjectionFallbackCacheUpToDate_ = true;
  return;
}
// ============================================================================
} // namespace Ginla
