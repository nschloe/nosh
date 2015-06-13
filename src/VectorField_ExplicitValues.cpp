// @HEADER
//
//    Query routines for a vector potential with explicitly given values.
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

#include "VectorField_ExplicitValues.hpp"

#include <map>
#include <string>

#include "Mesh.hpp"

namespace Nosh
{
namespace VectorField
{
// ============================================================================
ExplicitValues::
ExplicitValues(
    const Nosh::Mesh & mesh,
    const std::string & fieldName,
    const double mu
    ) :
  mu_(mu),
  edgeProjectionCache_(mesh.getMyEdges().size())
{
  // Initialize the cache.
  const std::vector<edge> edges = mesh.getMyEdges();

  const VectorFieldType & coordsField = mesh.getNodeField("coordinates");

  const VectorFieldType & dataField = mesh.getNodeField(fieldName);

  // Loop over all edges and create the cache.
  for (std::size_t k = 0; k < edges.size(); k++) {
    // Approximate the value at the midpoint of the edge
    // by the average of the values at the adjacent nodes.
    Eigen::Vector3d av = 0.5 * (
        mesh.getNodeValue(dataField, std::get<0>(edges[k]))
        + mesh.getNodeValue(dataField, std::get<1>(edges[k]))
        );

    // Extract the nodal coordinates.
    Eigen::Vector3d myEdge =
      mesh.getNodeValue(coordsField, std::get<1>(edges[k]))
      - mesh.getNodeValue(coordsField, std::get<0>(edges[k]));

    edgeProjectionCache_[k] = av.dot(myEdge);
  }

// TODO resurrect this
//#ifndef NDEBUG
#if 0
  // Do a quick sanity check for the edgeProjectionCache_.  It happens too
  // often that the reader elements aren't specified correctly and stk_io
  // *silently* "reads" only zeros.  Use the fake logical "isNonzeroLocal"
  // since Teuchos::Comm<int> doesn't have logical any() or all() operations.
  bool isZeroLocal = true;
  for (std::size_t k = 0; k < edgeProjectionCache_.size(); k++) {
    if (fabs(edgeProjectionCache_[k]) > 1.0e-10) {
      isZeroLocal = false;
      break;
    }
  }
  bool isZeroGlobal;
  Teuchos::reduceAll(
      *mesh.getComm(),
      Teuchos::REDUCE_AND,
      1,
      &isZeroLocal,
      &isZeroGlobal
      );

  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      isZeroGlobal,
      "Field \"" << fieldName << "\" seems empty. Was it read correctly?"
      );
#endif

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
setParameters(const std::map<std::string, double> & params)
{
  mu_ = params.at("mu");
  return;
}
// ============================================================================
const std::map<std::string, double>
ExplicitValues::
getParameters() const
{
  return {{"mu", mu_}};
}
// ============================================================================
double
ExplicitValues::
getEdgeProjection(const unsigned int edgeIndex) const
{
  return mu_ * edgeProjectionCache_[edgeIndex];
}
// ============================================================================
double
ExplicitValues::
getDEdgeProjectionDp(
    const unsigned int edgeIndex,
    const std::string & dParamName
    ) const
{
  if (dParamName.compare("mu") == 0)
    return edgeProjectionCache_[edgeIndex];
  else
    return 0.0;
}
// ============================================================================
} // namespace VectorField
} // namespace Nosh
