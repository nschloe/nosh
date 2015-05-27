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

#include "nosh/VectorField_ConstantCurl.hpp"

#include <map>
#include <string>

#include <Epetra_Vector.h>

#include "nosh/StkMesh.hpp"

namespace Nosh
{
namespace VectorField
{
// ============================================================================
ConstantCurl::
ConstantCurl(
    const Teuchos::RCP<Nosh::StkMesh> &mesh,
    const Teuchos::RCP<Eigen::Vector3d> &b,
    const Teuchos::RCP<Eigen::Vector3d> &u
    ) :
  mesh_(mesh),
  b_(b),
  u_(u),
  rotatedBCache_(*b),
  rotatedBCacheAngle_(0.0),
  dRotatedBDThetaCache_(*b_),
  rotateddBdThetaCacheAngle_(0.0),
  edgeCache_(),
  edgeCacheUptodate_(false),
  mu_(0.0),
  theta_(0.0)
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(!mesh_.is_null());
  TEUCHOS_ASSERT(!b_.is_null());
#endif
  TEUCHOS_TEST_FOR_EXCEPT_MSG(b_->dot(*b_) != 1.0,
                               "Curl vector not normalized: "
                               << "<b,b> = " << b->dot(*b) << "."
                               << std::endl
                            );
  if (!u_.is_null()) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(u_->dot(*u_) != 1.0,
                                 "Rotation vector not normalized: "
                                 << "<u,u> = " << u_->dot(*u_) << "."
                                 << std::endl
                              );
    this->dRotateDTheta_(dRotatedBDThetaCache_, *u_, 0.0);
  }

  edgeCache_ = Teuchos::ArrayRCP<Eigen::Vector3d>(mesh_->getEdgeNodes().size());

  return;
}
// ============================================================================
ConstantCurl::
~ConstantCurl()
{
}
// ============================================================================
void
ConstantCurl::
setParameters(const std::map<std::string, double> & params)
{
  mu_ = params.at("mu");
  theta_ = params.at("theta");
  return;
}
// ============================================================================
const std::map<std::string, double>
ConstantCurl::
getParameters() const
{
  return {{"mu", mu_}, {"theta", theta_}};
}
// ============================================================================
double
ConstantCurl::
getEdgeProjection(const unsigned int edgeIndex) const
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
  if (!edgeCacheUptodate_)
    this->initializeEdgeCache_();

  if (rotatedBCacheAngle_ != theta_) {
    rotatedBCache_ = *b_;
    this->rotate_(rotatedBCache_, *u_, theta_);
    rotatedBCacheAngle_ = theta_;
  }

  return mu_ * rotatedBCache_.dot(edgeCache_[edgeIndex]);
}
// ============================================================================
double
ConstantCurl::
getDEdgeProjectionDp(
    const unsigned int edgeIndex,
    const std::string & paramName
    ) const
{
  // Update caches.
  if (!edgeCacheUptodate_)
    this->initializeEdgeCache_();

  if (rotatedBCacheAngle_ != theta_) {
    rotatedBCache_ = *b_;
    this->rotate_(rotatedBCache_, *u_, theta_);
    rotatedBCacheAngle_ = theta_;
  }

  if (rotateddBdThetaCacheAngle_ != theta_) {
    dRotatedBDThetaCache_ = *b_;
    this->dRotateDTheta_(dRotatedBDThetaCache_, *u_, theta_);
    rotateddBdThetaCacheAngle_ = theta_;
  }

  if (paramName.compare("mu") == 0) {
    return rotatedBCache_.dot(edgeCache_[edgeIndex]);
  } else if (paramName.compare("theta") == 0) {
    return mu_ * dRotatedBDThetaCache_.dot(edgeCache_[edgeIndex]);
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "Illegal parameter \"" << paramName << "\"."
        );
  }
}
// ============================================================================
void
ConstantCurl::
rotate_(
    Eigen::Vector3d &v,
    const Eigen::Vector3d &u,
    const double theta
    ) const
{
  // Rotate a vector \c v by the angle \c theta in the plane perpendicular
  // to the axis given by \c u.
  // Refer to
  // http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
  double sinTheta = sin(theta);
  double cosTheta = cos(theta);

  if (sinTheta != 0.0) {
    Eigen::Vector3d vOld = v;

    // cos(theta) * I * v
    v *= cosTheta;

    // + sin(theta) u\cross v
    // Instead of what we have here,
    // we'd much rather write
    //   v += sinTheta * this->crossProduct_(u, vOld);
    // or do something like a DAXPY.
    // However, the Teuchos::SerialDenseVector doesn't have
    // that capability.
    Eigen::Vector3d tmp = this->crossProduct_(u, vOld);
    tmp *= sinTheta;
    v += tmp;

    // + (1-cos(theta)) (u*u^T) * v
    tmp = u;
    tmp *= (1.0-cosTheta) * u.dot(vOld);
    v += tmp;
  }

  return;
}
// ============================================================================
void
ConstantCurl::
dRotateDTheta_(
    Eigen::Vector3d &v,
    const Eigen::Vector3d &u,
    const double theta
    ) const
{
  // Incremental change of the rotation of a vector v around the axis u
  // by the angle theta.
  // Compare with method above.
  double sinTheta = sin(theta);
  double cosTheta = cos(theta);

  Eigen::Vector3d vOld = v;

  // -sin(theta) * I * v
  v *= -sinTheta;

  // + cos(theta) u\cross v
  //
  // Instead of what we have here,
  // we'd much rather write
  //   v += sinTheta * this->crossProduct_(u, vOld);
  // or do something like a DAXPY.
  // However, the Teuchos::SerialDenseVector doesn't have
  // that capability.
  Eigen::Vector3d tmp = this->crossProduct_(u, vOld);
  tmp *= cosTheta;
  v += tmp;

  // + (1+sin(theta)) (u*u^T) * v
  tmp = u;
  tmp *= (1.0+sinTheta) * u.dot(vOld);
  v += tmp;

  return;
}
// ============================================================================
Eigen::Vector3d
ConstantCurl::
crossProduct_(
    const Eigen::Vector3d u,
    const Eigen::Vector3d v
    ) const
{
  Eigen::Vector3d uXv(3);
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
  TEUCHOS_ASSERT(!mesh_.is_null());
#endif
  const Teuchos::Array<edge> edges = mesh_->getEdgeNodes();

  // Loop over all edges and create the cache.
  for (auto k = 0; k < edges.size(); k++) {
    const Eigen::Vector3d & node0Coords =
      mesh_->get3dVectorFieldNonconst(std::get<0>(edges[k]), "coordinates");
    const Eigen::Vector3d & node1Coords =
      mesh_->get3dVectorFieldNonconst(std::get<1>(edges[k]), "coordinates");

    // edgeMidpoint x edge = 0.5 (a+b) x (a-b) = b x a
    edgeCache_[k] = 0.5 * node0Coords.cross(node1Coords);
  }

  edgeCacheUptodate_ = true;

  return;
}
// ============================================================================
} // namespace VectorField
} // namespace Nosh
