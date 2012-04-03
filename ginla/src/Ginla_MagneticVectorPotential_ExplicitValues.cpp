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
ExplicitValues( const Teuchos::RCP<Ginla::StkMesh> &mesh,
                const Teuchos::RCP<const Epetra_MultiVector> &mvp,
                double mu
                ) :
  mesh_( mesh ),
  mvp_( mvp ),
  mu_( mu ),
  mvpEdgeMidpointProjectionCache_( Teuchos::ArrayRCP<double>() ),
  mvpEdgeMidpointProjectionCacheUptodate_( false )
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif

  mvpEdgeMidpointProjectionCache_ =
    Teuchos::ArrayRCP<double>( mesh_->getEdgeNodes().size() );

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
setParameters( const LOCA::ParameterVector &p )
{
  if (p.isParameter( "mu" ))
    mu_ = p.getValue( "mu" );

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
  TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
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

  const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > edges =
    mesh_->getEdgeNodes();

  // Loop over all edges and create the cache.
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
    // Get the two end points.
    Teuchos::Tuple<int,2> lid;
    lid[0] = mvp_->Map().LID( edges[k][0]->identifier() - 1 );
    lid[1] = mvp_->Map().LID( edges[k][1]->identifier() - 1 );

#ifdef _DEBUG_
    TEUCHOS_TEST_FOR_EXCEPT_MSG( lid[0] < 0,
                         "The global index " <<
                         edges[k][0]->identifier() - 1
                         << " does not seem to be present on this node." );
    TEUCHOS_TEST_FOR_EXCEPT_MSG( lid[1] < 0,
                         "The global index " <<
                         edges[k][1]->identifier() - 1
                         << " does not seem to be present on this node." );
#endif

    // Approximate the value at the midpoint of the edge
    // by the average of the values at the adjacent nodes.
    DoubleVector mvpEdgeMidpoint( 3 );
    for (int i=0; i<3; i++ )
      mvpEdgeMidpoint[i] = 0.5
                         * ((*(*mvp_)(i))[lid[0]] + (*(*mvp_)(i))[lid[1]]);

    // extract the nodal coordinates
    DoubleVector edge = mesh_->getNodeCoordinates( edges[k][1] );
    edge -= mesh_->getNodeCoordinates( edges[k][0] );

    mvpEdgeMidpointProjectionCache_[k] = edge.dot( mvpEdgeMidpoint );
  }

  mvpEdgeMidpointProjectionCacheUptodate_ = true;
  return;
}
// ============================================================================
} // namespace MagneticVectorPotential
} // namespace Ginla
