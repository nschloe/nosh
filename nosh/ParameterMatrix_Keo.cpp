// @HEADER
//
//    Builds the kinetic energy operator and its variants.
//    Copyright (C) 2010--2012  Nico Schl\"omer
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
// =============================================================================
// includes
#include "nosh/ParameterMatrix_Keo.hpp"

#include <map>
#include <string>

#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Virtual.hpp"
#include "nosh/VectorField_Virtual.hpp"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>

#include <ml_epetra_preconditioner.h>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif

namespace Nosh
{
namespace ParameterMatrix
{
// =============================================================================
Keo::
Keo(
    const std::shared_ptr<const Nosh::StkMesh> &mesh,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
    const std::shared_ptr<Nosh::VectorField::Virtual> &mvp
   ):
  Virtual(mesh),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  keoFillTime_(Teuchos::TimeMonitor::getNewTimer("Nosh: Keo::fill_")),
#endif
  thickness_(thickness),
  mvp_(mvp),
  alphaCache_(),
  alphaCacheUpToDate_(false)
{
}
// =============================================================================
Keo::
~Keo()
{
}
// =============================================================================
std::shared_ptr<Virtual>
Keo::
clone() const
{
  return std::shared_ptr<Keo>(new Keo(*this));
}
// =============================================================================
const std::map<std::string, double>
Keo::
getParameters() const
{
  return mvp_->getParameters();
}
// =============================================================================
double
Keo::
integrate1d_(
    const Nosh::VectorField::Virtual & f,
    const Eigen::Vector3d & x0,
    const Eigen::Vector3d & x1
    ) const
{
  if (f.degree() < 2) {
    return (x0 - x1).dot(mvp_->eval(0.5 * (x0 + x1)));
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "Cannot integrate functions of degree " << f.degree() << "."
        );
  }
}
// =============================================================================
void
Keo::
refill_(const std::map<std::string, double> & params)
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*keoFillTime_);
#endif

  mvp_->setParameters(params);

  // Zero-out the matrix.
  TEUCHOS_ASSERT_EQUALITY(0, this->PutScalar(0.0));

#ifndef NDEBUG
  TEUCHOS_ASSERT(mesh_);
  TEUCHOS_ASSERT(thickness_);
  TEUCHOS_ASSERT(mvp_);
#endif

  const std::vector<edge> edges = mesh_->getEdgeNodes();
  if (!alphaCacheUpToDate_) {
    this->buildAlphaCache_(edges, mesh_->getEdgeCoefficients());
  }

  double v[3];
  const VectorFieldType & coordsField = mesh_->getNodeField("coordinates");
  // Loop over all edges.
  for (auto k = 0; k < edges.size(); k++) {
    // Compute the integral
    //
    //    I = \int_{x0}^{xj} (xj-x0).A(x) / |xj-x0| dx
    //
    // numerically by the midpoint rule, i.e.,
    //
    //    I ~ |xj-x0| * (xj-x0) . A(0.5*(xj+x0)) / |xj-x0|.
    //
    // Project vector field onto the edge.
    // Instead of first computing the projection over the normalized edge
    // and then multiply it with the edge length, don't normalize the
    // edge vector.
    // We'd like to insert the 2x2 matrix
    //
    //     [   alpha                 , - alpha * exp(-IM * aInt) ]
    //     [ - alpha * exp(IM * aInt),   alpha                   ]
    //
    // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
    // that shares and edge.
    // Do that now, just blockwise for real and imaginary part.

    const double aInt = mvp_->getEdgeProjection(k);

    //// ----
    //// get edge coords (cache this)
    //// Now: midpoint rule
    //const Eigen::Vector3d x0 =
    //  mesh_->getNodeValue(coordsField, std::get<0>(edges[k]));
    //const Eigen::Vector3d x1 =
    //  mesh_->getNodeValue(coordsField, std::get<1>(edges[k]));
    //const double aInt2 = integrate1d_(*mvp_, x0, x1);

    //TEUCHOS_TEST_FOR_EXCEPT_MSG(
    //    fabs(aInt - aInt2) > 1.0e-10,
    //    aInt << " != " << aInt2 << " diff  = " << aInt - aInt2
    //   );

    // ----
    // Fill v with
    // Re(-exp(i Aint))
    // Im(-exp(i Aint))
    // 1.0
    double sinAInt, cosAInt;
    sincos(aInt, &sinAInt, &cosAInt);
    v[0] = -cosAInt * alphaCache_[k];
    v[1] = -sinAInt * alphaCache_[k];
    v[2] = alphaCache_[k];
    double ain [] = {
      v[2],  0.0,   v[0], v[1],
       0.0,  v[2], -v[1], v[0],
      v[0], -v[1],  v[2],  0.0,
      v[1],  v[0],   0.0, v[2]
    };
    int ierr = this->SumIntoGlobalValues(
        mesh_->globalIndexCache[k],
        Epetra_SerialDenseMatrix(View, ain, 4, 4 ,4)
        );
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(0, ierr);
#endif
  }

  // calls FillComplete by default
  int ierr = this->GlobalAssemble();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(0, ierr);
#endif
  return;
}
// =============================================================================
void
Keo::
buildAlphaCache_(
    const std::vector<edge> & edges,
    const std::vector<double> & edgeCoefficients
    ) const
{
  // This routine serves the one and only purpose of caching the
  // thickness average. The cache is used in every call to this->fill().
  // This is somewhat problematic since the map of V is principally
  // not known here. Also, it is typically a nonoverlapping map whereas
  // some edges do sit on a processor boundary, so actually the values
  // of V are needed in an overlapping map.
  // Fair enough. Let's distribute the vales of V to an overlapping
  // map here.
  alphaCache_ = std::vector<double>(edges.size());

  std::map<std::string, double> dummy;
  const Epetra_Vector thicknessValues = thickness_->getV(dummy);

  std::shared_ptr<const Epetra_Map> overlapMap = mesh_->getNodesOverlapMap();
  // We need to make sure that thicknessValues are distributed on
  // the overlap map.
  // Make sure to use Import here instead of Export as the vector
  // that we want to build is overlapping, "larger". If the "smaller",
  // non-overlapping vector is exported, only the values on the overlap
  // would only be set on one processor.
  Epetra_Vector thicknessOverlap(*overlapMap);
  Epetra_Import importer(*overlapMap, thicknessValues.Map());
  TEUCHOS_ASSERT_EQUALITY(
      0,
      thicknessOverlap.Import(thicknessValues, importer, Insert)
      );

  int gid0, gid1;
  int lid0, lid1;
  for (auto k = 0; k < edges.size(); k++) {
    // Get the ID of the edge endpoints in the map of
    // getV(). Well...
    gid0 = mesh_->gid(std::get<0>(edges[k]));
    lid0 = overlapMap->LID(gid0);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid0 < 0,
        "The global index " << gid0
        << " does not seem to be present on this node."
       );
#endif
    gid1 = mesh_->gid(std::get<1>(edges[k]));
    lid1 = overlapMap->LID(gid1);
#ifndef NDEBUG
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lid1 < 0,
        "The global index " << gid1
        << " does not seem to be present on this node."
       );
#endif
    // Update cache.
    alphaCache_[k] = edgeCoefficients[k]
      * 0.5 * (thicknessOverlap[lid0] + thicknessOverlap[lid1]);
  }

  alphaCacheUpToDate_ = true;

  return;
}
// =============================================================================
}  // namespace ParameterMatrix
}  // namespace Nosh
