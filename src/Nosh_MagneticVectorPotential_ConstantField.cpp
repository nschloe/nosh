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

#include "Nosh_MagneticVectorPotential_ConstantField.hpp"
#include "Nosh_StkMesh.hpp"

#include <Epetra_Vector.h>

namespace Nosh {
namespace MagneticVectorPotential {
// ============================================================================
ConstantField::
ConstantField( const Teuchos::RCP<Nosh::StkMesh> &mesh,
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
  edgeCacheUptodate_( false )
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !b_.is_null() );
#endif
  TEUCHOS_TEST_FOR_EXCEPT_MSG( b_->dot( *b_ ) != 1.0,
                       "Magnetic field vector not normalized: "
                       << "<b,b> = " << b->dot( *b ) << "."
                       << std::endl
                       );
  if ( !u_.is_null() )
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG( u_->dot( *u_ ) != 1.0,
                         "Rotation vector not normalized: "
                         << "<u,u> = " << u_->dot( *u_ ) << "."
                         << std::endl
                         );
    // If compiled with GNU (and maybe other compilers), we could use
    // sincos() here to comute sin and cos simultaneously.
    // PGI, for one, doesn't support sincos, though.
    double sinTheta = sin(theta_);
    double cosTheta = cos(theta_);
    //sincos( theta_, &sinTheta, &cosTheta );
    rotatedB_ = this->rotate_( *b_, *u_, sinTheta, cosTheta );
    dRotatedBDTheta_ = this->dRotateDTheta_( *b_, *u_, sinTheta, cosTheta );
  }


  edgeCache_ = Teuchos::ArrayRCP<DoubleVector>(mesh_->getEdgeNodes().size());

  return;
}
// ============================================================================
ConstantField::
~ConstantField()
{
}
// ============================================================================
DoubleVector
ConstantField::
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
ConstantField::
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
ConstantField::
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
ConstantField::
setParameters( const LOCA::ParameterVector &p )
{
  if (p.isParameter( "mu" ))
    mu_ = p.getValue( "mu" );

  if (p.isParameter( "theta" ))
  {
    theta_ = p.getValue( "theta" );
    // If compiled with GNU (and maybe other compilers), we could use
    // sincos() here to comute sin and cos simultaneously.
    // PGI, for one, doesn't support sincos,
    double sinTheta = sin(theta_);
    double cosTheta = cos(theta_);
    //sincos( theta_, &sinTheta, &cosTheta );
    rotatedB_ = this->rotate_( *b_, *u_, sinTheta, cosTheta );
    dRotatedBDTheta_ = this->dRotateDTheta_( *b_, *u_, sinTheta, cosTheta );
  }

  return;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
ConstantField::
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
ConstantField::
getAEdgeMidpointProjection(const unsigned int edgeIndex
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
ConstantField::
getdAdPEdgeMidpointProjection(const unsigned int edgeIndex,
                              const unsigned int parameterIndex
                              ) const
{
  if ( !edgeCacheUptodate_ )
    this->initializeEdgeCache_();

  if (parameterIndex == 0) // mu
    return rotatedB_.dot( edgeCache_[edgeIndex] );
  else if (parameterIndex == 1) // theta
    return mu_ * dRotatedBDTheta_.dot( edgeCache_[edgeIndex] );
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                                "Unknown parameter index " << parameterIndex << ".");
}
// ============================================================================
void
ConstantField::
initializeEdgeCache_() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > edges =
    mesh_->getEdgeNodes();

  // Loop over all edges and create the cache.
  for (unsigned int k=0; k<edges.size(); k++)
  {
    const DoubleVector & node0Coords =
      mesh_->getNodeCoordinatesNonconst( edges[k][0] );
    const DoubleVector & node1Coords =
      mesh_->getNodeCoordinatesNonconst( edges[k][1] );

    // edgeMidpoint x edge = 0.5 (a+b) x (a-b) = b x a
    edgeCache_[k] = this->crossProduct_( node0Coords, node1Coords );
    TEUCHOS_ASSERT_EQUALITY( 0, edgeCache_[k].scale( 0.5 ) );
  }

  edgeCacheUptodate_ = true;

  return;
}
// ============================================================================
} // namespace MagneticVectorPotential
} // namespace Nosh
