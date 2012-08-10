// @HEADER
//
//    Query routines for a vector potential with explicitly given values.
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

#include "Nosh_VectorField_ExplicitValues.hpp"
#include "Nosh_StkMesh.hpp"

#include <Epetra_Vector.h>

namespace Nosh {
namespace VectorField {
// ============================================================================
ExplicitValues::
ExplicitValues(const Teuchos::RCP<Nosh::StkMesh> &mesh,
               const Teuchos::RCP<const Epetra_MultiVector> &values,
               const double initMu
               ) :
  mesh_( mesh ),
  values_( values ),
  initMu_( initMu ),
  edgeProjectionCache_( Teuchos::ArrayRCP<double>() ),
  edgeProjectionCacheUptodate_( false )
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif

  edgeProjectionCache_ =
    Teuchos::ArrayRCP<double>( mesh_->getEdgeNodes().size() );

  return;
}
// ============================================================================
ExplicitValues::
~ExplicitValues()
{
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<double> >
ExplicitValues::
get_p_init() const
{
  Teuchos::RCP<Teuchos::Array<double> > p_init =
    Teuchos::rcp(new Teuchos::Array<double>(1));
  (*p_init)[0] = initMu_;
  return p_init;
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
ExplicitValues::
get_p_names() const
{
  Teuchos::RCP<Teuchos::Array<std::string> > p_names =
    Teuchos::rcp(new Teuchos::Array<std::string>(1));
  (*p_names)[0] = "mu";
  return p_names;
}
// ============================================================================
double
ExplicitValues::
getEdgeProjection(const unsigned int edgeIndex,
                  const Teuchos::Array<double> & params
                  ) const
{
  if ( !edgeProjectionCacheUptodate_ )
    this->initializeEdgeMidpointCache_();

  return params[0] * edgeProjectionCache_[edgeIndex];
}
// ============================================================================
double
ExplicitValues::
getDEdgeProjectionDp(const unsigned int edgeIndex,
                     const Teuchos::Array<double> & params,
                     const unsigned int parameterIndex
                     ) const
{
  TEUCHOS_ASSERT_EQUALITY(parameterIndex, 0);

  if ( !edgeProjectionCacheUptodate_ )
    this->initializeEdgeMidpointCache_();

  return edgeProjectionCache_[edgeIndex];
}
// ============================================================================
void
ExplicitValues::
initializeEdgeMidpointCache_() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !values_.is_null() );
#endif

  const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > edges =
    mesh_->getEdgeNodes();

  // Loop over all edges and create the cache.
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
    // Get the two end points.
    Teuchos::Tuple<int,2> lid;
    lid[0] = values_->Map().LID( edges[k][0]->identifier() - 1 );
    lid[1] = values_->Map().LID( edges[k][1]->identifier() - 1 );

#ifndef NDEBUG
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
    DoubleVector valueEdgeMidpoint( 3 );
    for (int i=0; i<3; i++ )
      valueEdgeMidpoint[i] = 0.5
                           * ((*(*values_)(i))[lid[0]] + (*(*values_)(i))[lid[1]]);

    // extract the nodal coordinates
    DoubleVector edge = mesh_->getNodeCoordinatesNonconst( edges[k][1] );
    edge -= mesh_->getNodeCoordinatesNonconst( edges[k][0] );

    //edgeProjectionCache_[k] = edge.dot( valueEdgeMidpoint );
    edgeProjectionCache_[k] = valueEdgeMidpoint.dot(edge);
  }

  edgeProjectionCacheUptodate_ = true;
  return;
}
// ============================================================================
} // namespace VectorField
} // namespace Nosh
