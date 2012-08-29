// @HEADER
//
//    Query routines for the vector potential associated with a constant curl field.
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

#include "Nosh_VectorField_ConstantCurl.hpp"
#include "Nosh_StkMesh.hpp"

#include <Epetra_Vector.h>

namespace Nosh {
namespace VectorField {
// ============================================================================
ConstantCurl::
ConstantCurl(const Teuchos::RCP<Nosh::StkMesh> &mesh,
             const Teuchos::RCP<DoubleVector> &b,
             const Teuchos::RCP<DoubleVector> &u
             ) :
  mesh_( mesh ),
  b_( b ),
  u_( u ),
  rotatedBCache_( *b ),
  rotatedBCacheAngle_(0.0),
  dRotatedBDThetaCache_(*b_),
  rotateddBdThetaCacheAngle_(0.0),
  edgeCache_(),
  edgeCacheUptodate_( false )
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !b_.is_null() );
#endif
  TEUCHOS_TEST_FOR_EXCEPT_MSG( b_->dot( *b_ ) != 1.0,
                       "Curl vector not normalized: "
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
    this->dRotateDTheta_(dRotatedBDThetaCache_, *u_, 0.0);
  }

  edgeCache_ = Teuchos::ArrayRCP<DoubleVector>(mesh_->getEdgeNodes().size());

  return;
}
// ============================================================================
ConstantCurl::
~ConstantCurl()
{
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
ConstantCurl::
get_p_names() const
{
  Teuchos::RCP<Teuchos::Array<std::string> > p_names =
    Teuchos::rcp(new Teuchos::Array<std::string>(1));
  (*p_names)[0] = "mu";
  (*p_names)[1] = "theta";
  return p_names;
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<double> >
ConstantCurl::
get_p_init() const
{
  Teuchos::RCP<Teuchos::Array<double> > p_init =
    Teuchos::rcp(new Teuchos::Array<double>(2));
  (*p_init)[0] = 0.0;
  (*p_init)[1] = 0.0;
  return p_init;
}
// ============================================================================
double
ConstantCurl::
getEdgeProjection(const unsigned int edgeIndex,
                  const Teuchos::Array<double> &params
                  ) const
{
  // A vector potential associated with the constant curl field RB is
  //
  //     A(X) = 0.5 * (RB x X).
  //
  // Projecting it onto an edge e yields
  //     e * A(X) = e * (0.5*(RB x X))
  //              = RB * (0.5*(X x e)).
  //
  // Evaluating this expression at the edge midpoint Em of the edge e gives
  //     e * A(Em) = RB * (0.5*(Em x e)).
  //
  // The vector 0.5*(Em x e) can be cached for the mesh.
  // This saves caching e and the edge midpoint separately and also avoids
  // computing the cross-products more than once if B changes.

  // Update caches.
  if ( !edgeCacheUptodate_ )
    this->initializeEdgeCache_();

  if (rotatedBCacheAngle_ != params[1])
  {
    rotatedBCache_ = *b_;
    this->rotate_(rotatedBCache_, *u_, params[1]);
    rotatedBCacheAngle_ = params[1];
  }

  return params[0] * rotatedBCache_.dot( edgeCache_[edgeIndex] );
}
// ============================================================================
double
ConstantCurl::
getDEdgeProjectionDp(const unsigned int edgeIndex,
                     const Teuchos::Array<double> &params,
                     const unsigned int parameterIndex
                     ) const
{
  // Update caches.
  if ( !edgeCacheUptodate_ )
    this->initializeEdgeCache_();

  if (rotatedBCacheAngle_ != params[1])
  {
    rotatedBCache_ = *b_;
    this->rotate_(rotatedBCache_, *u_, params[1]);
    rotatedBCacheAngle_ = params[1];
  }

  if (rotateddBdThetaCacheAngle_ != params[1])
  {
    dRotatedBDThetaCache_ = *b_;
    this->dRotateDTheta_(dRotatedBDThetaCache_, *u_, params[1]);
    rotateddBdThetaCacheAngle_ = params[1];
  }

  if (parameterIndex == 0) // mu
    return rotatedBCache_.dot( edgeCache_[edgeIndex] );
  else if (parameterIndex == 1) // theta
    return params[0] * dRotatedBDThetaCache_.dot( edgeCache_[edgeIndex] );
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                                "Unknown parameter index " << parameterIndex << ".");
}
// ============================================================================
void
ConstantCurl::
rotate_(DoubleVector &v,
        const DoubleVector &u,
        const double theta
        ) const
{
  // Rotate a vector \c v by the angle \c theta in the plane perpendicular
  // to the axis given by \c u.
  // Refer to
  // http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
  double sinTheta = sin(theta);
  double cosTheta = cos(theta);

  if ( sinTheta != 0.0 )
  {
    DoubleVector vOld = v;

    // cos(theta) * I * v
    v *= cosTheta;

    // + sin(theta) u\cross v
    // Instead of what we have here,
    // we'd much rather write
    //   v += sinTheta * this->crossProduct_( u, vOld );
    // or do something like a DAXPY.
    // However, the Teuchos::SerialDenseVector doesn't have
    // that capability.
    DoubleVector tmp = this->crossProduct_( u, vOld );
    tmp *= sinTheta;
    v += tmp;

    // + (1-cos(theta)) (u*u^T) * v
    tmp = u;
    tmp *= (1.0-cosTheta) * u.dot( vOld );
    v += tmp;
  }

  return;
}
// ============================================================================
void
ConstantCurl::
dRotateDTheta_(DoubleVector &v,
               const DoubleVector &u,
               const double theta
               ) const
{
  // Incremental change of the rotation of a vector v around the axis u
  // by the angle theta.
  // Compare with method above.
  double sinTheta = sin(theta);
  double cosTheta = cos(theta);

  DoubleVector vOld = v;

  // -sin(theta) * I * v
  v *= -sinTheta;

  // + cos(theta) u\cross v
  //
  // Instead of what we have here,
  // we'd much rather write
  //   v += sinTheta * this->crossProduct_( u, vOld );
  // or do something like a DAXPY.
  // However, the Teuchos::SerialDenseVector doesn't have
  // that capability.
  DoubleVector tmp = this->crossProduct_( u, vOld );
  tmp *= cosTheta;
  v += tmp;

  // + (1+sin(theta)) (u*u^T) * v
  tmp = u;
  tmp *= (1.0+sinTheta) * u.dot( vOld );
  v += tmp;

  return;
}
// ============================================================================
DoubleVector
ConstantCurl::
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
ConstantCurl::
initializeEdgeCache_() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > edges =
    mesh_->getEdgeNodes();

  // Loop over all edges and create the cache.
  for (unsigned int k=0; k<edges.size(); k++)
  {
    const DoubleVector & node0Coords =
      mesh_->getVectorFieldNonconst(edges[k][0],
                                    "coordinates", 3);
    const DoubleVector & node1Coords =
      mesh_->getVectorFieldNonconst(edges[k][1],
                                    "coordiates", 3);

    // edgeMidpoint x edge = 0.5 (a+b) x (a-b) = b x a
    edgeCache_[k] = this->crossProduct_( node0Coords, node1Coords );
    TEUCHOS_ASSERT_EQUALITY( 0, edgeCache_[k].scale( 0.5 ) );
  }

  edgeCacheUptodate_ = true;

  return;
}
// ============================================================================
} // namespace VectorField
} // namespace Nosh
