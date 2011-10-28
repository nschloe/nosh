/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2011  Nico Schl\"omer

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

#include "Ginla_MagneticVectorPotential.hpp"
#include "Ginla_StkMesh.hpp"

#include <Epetra_Vector.h>

namespace Ginla {
// ============================================================================
MagneticVectorPotential::
MagneticVectorPotential( const Teuchos::RCP<Ginla::StkMesh> & mesh,
                         const Teuchos::RCP<const Epetra_MultiVector>  & mvp,
                         double mu
                       ):
  mesh_( mesh ),
  mvp_( mvp ),
  mu_( mu ),
  edgeMidpointProjectionCache_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( mesh->getOwnedCells().size() ) ),
  edgeMidpointProjectionCacheUpToDate_( false )
{
}
// ============================================================================
MagneticVectorPotential::
~MagneticVectorPotential()
{
}
// ============================================================================
bool
MagneticVectorPotential::
setParameters( const LOCA::ParameterVector & p )
{
    bool valuesChanged = false;

    if (p.isParameter( "mu" ))
        if ( mu_ != p.getValue ( "mu" ) )
        {
            mu_ = p.getValue ( "mu" );
            valuesChanged = true;
        }

    return valuesChanged;
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
Teuchos::RCP<Point>
MagneticVectorPotential::
getA(const Point & x) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                               "getA(x) for general 'x' cannot be implemented as "
                               << "no general MVP info is provided to this class."
                             );
  Teuchos::RCP<Point> a;
  return a;
}
// ============================================================================
double
MagneticVectorPotential::
getAEdgeMidpointProjection( const unsigned int cellIndex,
                            const unsigned int edgeIndex
                          ) const
{
  if ( !edgeMidpointProjectionCacheUpToDate_ )
      this->initializeEdgeMidpointProjectionCache_();

  return mu_ * edgeMidpointProjectionCache_[cellIndex][edgeIndex];
}
// ============================================================================
void
MagneticVectorPotential::
initializeEdgeMidpointProjectionCache_() const
{
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // Loop over all edges and create the cache.
  // To this end, loop over all cells and the edges within the cell.
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation rel = (*cells[k]).relations();

      unsigned int numLocalNodes = rel.size();
      unsigned int cellDimension = mesh_->getCellDimension( numLocalNodes );
      // extract the nodal coordinates
      Teuchos::Array<Point> localNodes = mesh_->getNodeCoordinates( rel );

      edgeMidpointProjectionCache_[k] = Teuchos::ArrayRCP<double>( mesh_->getNumEdgesPerCell( cellDimension ) );

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      unsigned int edgeIndex = 0;
      Teuchos::Tuple<int,2> gid, lid;
      for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
      {
          const Point & node0 = localNodes[e0];
          gid[0] = (*rel[e0].entity()).identifier() - 1;
          lid[0] = mvp_->Map().LID( gid[0] );
          TEUCHOS_TEST_FOR_EXCEPT_MSG( lid[0] < 0,
                                       "The global index " << gid[0]
                                       << " does not seem to be present on this node." );
          for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
          {
              const Point & node1 = localNodes[e1];
              gid[1] = (*rel[e1].entity()).identifier() - 1;
              lid[1] = mvp_->Map().LID( gid[1] );
              TEUCHOS_TEST_FOR_EXCEPT_MSG( lid[1] < 0,
                                           "The global index " << gid[1]
                                           << " does not seem to be present on this node." );

              // Approximate the value at the midpoint of the edge
              // by the average of the values at the adjacent nodes.
              Point a;
              for (int i=0; i<a.size(); i++ )
                  a[i] = 0.5 * ( (*(*mvp_)(i))[lid[0]] + (*(*mvp_)(i))[lid[1]] );

              // Instead of first computing the projection over the normalized edge
              // and then multiply it with the edge length, don't normalize the
              // edge vector.
              edgeMidpointProjectionCache_[k][edgeIndex] = 0.0;
              for (int i=0; i<a.size(); i++ )
                  edgeMidpointProjectionCache_[k][edgeIndex] += ( node1[i] - node0[i] ) * a[i];

              edgeIndex++;
          }
      }
  }

  edgeMidpointProjectionCacheUpToDate_ = true;
  return;
}
// ============================================================================
} // namespace Ginla
