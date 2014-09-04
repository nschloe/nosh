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

#include "nosh/VectorField_ExplicitValues.hpp"

#include <map>
#include <string>

#include <Epetra_Vector.h>

#include "nosh/StkMesh.hpp"

namespace Nosh
{
namespace VectorField
{
// ============================================================================
ExplicitValues::
ExplicitValues(const Nosh::StkMesh & mesh,
               const std::string & fieldName,
               const double initMu
             ) :
  initMu_(initMu),
  edgeProjectionCache_(mesh.getEdgeNodes().size())
{
  typedef std::tuple<stk::mesh::Entity,stk::mesh::Entity> edge;
  // Initialize the cache.
  const Teuchos::Array<edge> edges = mesh.getEdgeNodes();

  // Loop over all edges and create the cache.
  for (Teuchos::Array<edge>::size_type k = 0;
       k < edges.size();
       k++) {
    // Approximate the value at the midpoint of the edge
    // by the average of the values at the adjacent nodes.
    DoubleVector av = mesh.getVectorFieldNonconst(
        std::get<0>(edges[k]),
        fieldName, 3
        );
    av += mesh.getVectorFieldNonconst(
        std::get<1>(edges[k]),
        fieldName, 3
        );
    av *= 0.5;

    // Extract the nodal coordinates.
    DoubleVector myEdge =
      mesh.getVectorFieldNonconst(
          std::get<1>(edges[k]),
          "coordinates", 3
          );
    myEdge -=
      mesh.getVectorFieldNonconst(
          std::get<0>(edges[k]),
          "coordinates", 3
          );

    edgeProjectionCache_[k] = av.dot(myEdge);
  }

  // Do a quick sanity check for the edgeProjectionCache_.
  // It happens too effing often that the reader elements aren't specified
  // correctly and stk_io *silently* "reads" only zeros.
  // Use the fake logical "isNonzeroLocal" since Epetra_Comm doesn't have
  // logical any() or all() operations.
  int isNonzeroLocal = 0;
  for (int k = 0; k < edges.size(); k++) {
    if (fabs(edgeProjectionCache_[k]) > 1.0e-10) {
      isNonzeroLocal = 1;
      break;
    }
  }
  int isNonzeroGlobal;
  mesh.getComm().SumAll(&isNonzeroLocal, &isNonzeroGlobal, 1);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(isNonzeroGlobal == 0,
                              "Field \"" << fieldName << "\" seems empty. "
                              << "Was it read correctly?");

  return;
}
// ============================================================================
ExplicitValues::
~ExplicitValues()
{
}
// ============================================================================
const std::map<std::string, double>
ExplicitValues::
getInitialParameters() const
{
  std::map<std::string, double> m;
  m["mu"] = initMu_;
  return m;
}
// ============================================================================
double
ExplicitValues::
getEdgeProjection(const unsigned int edgeIndex,
                  const std::map<std::string, double> & params
                ) const
{
  std::map<std::string, double>::const_iterator it = params.find("mu");
  TEUCHOS_ASSERT(it != params.end());
  return it->second * edgeProjectionCache_[edgeIndex];
}
// ============================================================================
double
ExplicitValues::
getDEdgeProjectionDp(const unsigned int edgeIndex,
                     const std::map<std::string, double> & params,
                     const std::string & dParamName
                   ) const
{
  (void) params;
  if (dParamName.compare("mu") == 0)
    return edgeProjectionCache_[edgeIndex];
  else
    return 0.0;
}
// ============================================================================
} // namespace VectorField
} // namespace Nosh
