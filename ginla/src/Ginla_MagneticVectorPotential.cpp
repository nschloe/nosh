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
                         double thetaX,
                         double thetaY,
                         double thetaZ
                       ):
  mesh_( mesh ),
  mvp_( mvp ),
  mu_( mu ),
  thetaX_( thetaX ),
  thetaY_( thetaY ),
  thetaZ_( thetaZ ),
  mvpEdgeMidpoint_( Teuchos::ArrayRCP<DoubleVector>() ),
  edges_( Teuchos::ArrayRCP<DoubleVector>() ),
  mvpEdgeMidpointUpToDate_( false ),
  mvpEdgeMidpointFallback_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> >() ),
  edgesFallback_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> >() ),
  mvpEdgeMidpointFallbackUpToDate_( false )
{
    TEUCHOS_ASSERT( !mesh_.is_null() );
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

    if (p.isParameter( "thetaX" ))
            thetaX_ = p.getValue ( "thetaX" );

    if (p.isParameter( "thetaY" ))
            thetaY_ = p.getValue ( "thetaY" );

    if (p.isParameter( "thetaZ" ))
            thetaZ_ = p.getValue ( "thetaZ" );

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
  p->addParameter( "thetaX", thetaX_ );
  p->addParameter( "thetaY", thetaY_ );
  p->addParameter( "thetaZ", thetaZ_ );

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

  return mu_ * edges_[edgeIndex].dot( mvpEdgeMidpoint_[edgeIndex] );
}
// ============================================================================
double
MagneticVectorPotential::
getdAdMuEdgeMidpointProjection( const unsigned int edgeIndex
                              ) const
{
  if ( !mvpEdgeMidpointUpToDate_ )
      this->initializeMvpEdgeMidpointCache_();

  return edges_[edgeIndex].dot( mvpEdgeMidpoint_[edgeIndex] );
}
// ============================================================================
void
MagneticVectorPotential::
initializeMvpEdgeMidpointCache_() const
{
  TEUCHOS_ASSERT( !mesh_.is_null() );
  std::vector<stk::mesh::Entity*> edges = mesh_->getOverlapEdges();

  TEUCHOS_ASSERT( !mvp_.is_null() );
  // Loop over all edges and create the cache.
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
      // Get the two end points.
      const stk::mesh::PairIterRelation endPoints =
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
      mvpEdgeMidpoint_[k] = Teuchos::SerialDenseVector<int,double>(3);
      for (int i=0; i<3; i++ )
          mvpEdgeMidpoint_[k][i] = 0.5 * ( (*(*mvp_)(i))[lid[0]]
                                         + (*(*mvp_)(i))[lid[1]] );

      // extract the nodal coordinates
      Teuchos::ArrayRCP<Teuchos::SerialDenseVector<int,double> > localNodes =
          mesh_->getNodeCoordinates( endPoints );
      TEUCHOS_ASSERT_EQUALITY( localNodes.size(), 2 );
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

  return mu_ * edgesFallback_[cellIndex][edgeIndex].dot( mvpEdgeMidpointFallback_[cellIndex][edgeIndex] );
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

  return edgesFallback_[cellIndex][edgeIndex].dot( mvpEdgeMidpointFallback_[cellIndex][edgeIndex] );
}
// ============================================================================
void
MagneticVectorPotential::
initializeMvpEdgeMidpointFallback_() const
{
  TEUCHOS_ASSERT( !mesh_.is_null() );
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // Loop over all edges and create the cache.
  // To this end, loop over all cells and the edges within the cell.
  TEUCHOS_ASSERT( !mvp_.is_null() );
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation localNodes = (*cells[k]).relations( mesh_->getMetaData()->node_rank() );

      unsigned int numLocalNodes = localNodes.size();
      unsigned int cellDimension = mesh_->getCellDimension( numLocalNodes );
      // extract the nodal coordinates
      Teuchos::ArrayRCP<Teuchos::SerialDenseVector<int,double> > localNodeCoords =
          mesh_->getNodeCoordinates( localNodes );

      mvpEdgeMidpointFallback_[k] =
          Teuchos::ArrayRCP<DoubleVector>( mesh_->getNumEdgesPerCell( cellDimension ) );
      edgesFallback_[k] =
          Teuchos::ArrayRCP<DoubleVector>( mesh_->getNumEdgesPerCell( cellDimension ) );

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      unsigned int edgeIndex = 0;
      Teuchos::Tuple<int,2> gid, lid;
      for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
      {
          gid[0] = (*localNodes[e0].entity()).identifier() - 1;
          lid[0] = mvp_->Map().LID( gid[0] );
          TEST_FOR_EXCEPT_MSG( lid[0] < 0,
                               "The global index " << gid[0]
                               << " does not seem to be present on this node." );
          for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
          {
              gid[1] = (*localNodes[e1].entity()).identifier() - 1;
              lid[1] = mvp_->Map().LID( gid[1] );
              TEST_FOR_EXCEPT_MSG( lid[1] < 0,
                                   "The global index " << gid[1]
                                   << " does not seem to be present on this node." );

              // Approximate the value at the midpoint of the edge
              // by the average of the values at the adjacent nodes.
              mvpEdgeMidpointFallback_[k][edgeIndex] = Teuchos::SerialDenseVector<int,double>(3);
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
