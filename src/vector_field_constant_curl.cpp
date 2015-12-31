// @HEADER
//
//    Query routines for the vector potential associated with a constant curl field.
//    Copyright (C) 2012  Nico Schlömer
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

#include "vector_field_constant_curl.hpp"

#include <map>
#include <string>

#include "mesh.hpp"

namespace nosh
{
namespace vector_field
{
// ============================================================================
constantCurl::
constantCurl(
    const std::shared_ptr<nosh::mesh> &mesh,
    const std::shared_ptr<Eigen::Vector3d> &b,
    const std::shared_ptr<Eigen::Vector3d> &u
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
  TEUCHOS_ASSERT(mesh_);
  TEUCHOS_ASSERT(b_);
#endif
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      b_->dot(*b_) != 1.0,
      "Curl vector not normalized: <b,b> = " << b->dot(*b) << "." << std::endl
      );
  if (u_) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        u_->dot(*u_) != 1.0,
        "Rotation vector not normalized: <u,u> = " << u_->dot(*u_) << "."
        << std::endl
        );
    this->dRotateDTheta_(dRotatedBDThetaCache_, *u_, 0.0);
  }

  edgeCache_ = std::vector<Eigen::Vector3d>(mesh_->my_edges().size());

  return;
}
// ============================================================================
constantCurl::
~constantCurl()
{
}
// ============================================================================
void
constantCurl::
set_parameters(const std::map<std::string, double> & params)
{
  mu_ = params.at("mu");
  theta_ = params.at("theta");
  return;
}
// ============================================================================
const std::map<std::string, double>
constantCurl::
get_parameters() const
{
  return {{"mu", mu_}, {"theta", theta_}};
}
// ============================================================================
double
constantCurl::
get_edge_projection(const unsigned int edge_index) const
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
  if (!edgeCacheUptodate_) {
    this->initializeEdgeCache_();
  }

  if (rotatedBCacheAngle_ != theta_) {
    rotatedBCache_ = *b_;
    this->rotate_(rotatedBCache_, *u_, theta_);
    rotatedBCacheAngle_ = theta_;
  }

  return mu_ * rotatedBCache_.dot(edgeCache_[edge_index]);
}
// ============================================================================
double
constantCurl::
get_d_edge_projection_dp(
    const unsigned int edge_index,
    const std::string & param_name
    ) const
{
  // Update caches.
  if (!edgeCacheUptodate_) {
    this->initializeEdgeCache_();
  }

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

  if (param_name.compare("mu") == 0) {
    return rotatedBCache_.dot(edgeCache_[edge_index]);
  } else if (param_name.compare("theta") == 0) {
    return mu_ * dRotatedBDThetaCache_.dot(edgeCache_[edge_index]);
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "Illegal parameter \"" << param_name << "\"."
        );
  }
}
// ============================================================================
void
constantCurl::
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
  double sinTheta;
  double cosTheta;
  sincos(theta, &sinTheta, &cosTheta);

  if (sinTheta != 0.0) {
    Eigen::Vector3d vOld = v;

    // cos(theta) * I * v
    // + sin(theta) u\cross v
    // + (1-cos(theta)) (u*u^T) * v
    v = cosTheta * vOld
      + sinTheta * u.cross(vOld)
      + (1.0-cosTheta) * u.dot(vOld) * u;
  }

  return;
}
// ============================================================================
void
constantCurl::
dRotateDTheta_(
    Eigen::Vector3d &v,
    const Eigen::Vector3d &u,
    const double theta
    ) const
{
  // Incremental change of the rotation of a vector v around the axis u
  // by the angle theta.
  // Compare with method above.
  double sinTheta;
  double cosTheta;
  sincos(theta, &sinTheta, &cosTheta);

  Eigen::Vector3d vOld = v;

  // -sin(theta) * I * v
  // + cos(theta) u\cross v
  // + (1+sin(theta)) (u*u^T) * v
  v = -sinTheta * vOld
    + cosTheta * u.cross(vOld)
    + (1.0+sinTheta) * u.dot(vOld) * u;

  return;
}
// ============================================================================
void
constantCurl::
initializeEdgeCache_() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(mesh_);
#endif
  const std::vector<edge> edges = mesh_->my_edges();

#if 0
  const vector_fieldType & coords_field = mesh_->get_node_field("coordinates");

  // Loop over all edges and create the cache.
  for (std::size_t k = 0; k < edges.size(); k++) {
    const Eigen::Vector3d & node0Coords =
      mesh_->get_node_value(coords_field, std::get<0>(edges[k]));
    const Eigen::Vector3d & node1Coords =
      mesh_->get_node_value(coords_field, std::get<1>(edges[k]));

    // edge_midpoint x edge = 0.5 (a+b) x (a-b) = b x a
    edgeCache_[k] = 0.5 * node0Coords.cross(node1Coords);
  }

  edgeCacheUptodate_ = true;
#endif

  return;
}
// ============================================================================
} // namespace vector_field
} // namespace nosh
