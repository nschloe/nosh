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

#include "Nosh_ScalarField_ExplicitValues.hpp"
#include "Nosh_StkMesh.hpp"

#include <Epetra_Vector.h>

namespace Nosh {
namespace ScalarField {
// ============================================================================
ExplicitValues::
ExplicitValues(const Nosh::StkMesh & mesh,
               const std::string & fieldName
               ) :
  nodeValues_(mesh.getOwnedNodes().size())
{
  // Build the nodeValues_ vector for the field so the values can later
  // be fetched by their integer nodeIndex into the field.
  const std::vector<stk::mesh::Entity*> & ownedNodes =
    mesh.getOwnedNodes();

  for (int k=0; k<ownedNodes.size(); k++)
    nodeValues_[k] = mesh.getScalarFieldNonconst(ownedNodes[k], fieldName);

  return;
}
// ============================================================================
ExplicitValues::
~ExplicitValues()
{
}
// ============================================================================
const std::map<std::string,double>
ExplicitValues::
getParameters() const
{
  std::map<std::string,double> m;
  m["beta"] = 1.0;
  return m;
}
// ============================================================================
double
ExplicitValues::
getV(const unsigned int nodeIndex,
     const std::map<std::string,double> & params
     ) const
{
  double val = nodeValues_[nodeIndex];

  // Find the value of "beta" and use it as a factor.
  std::map<std::string, double>::const_iterator it = params.find("beta");
  TEUCHOS_ASSERT(it != params.end());
  val *= it->second;

  return val;
}
// ============================================================================
double
ExplicitValues::
getdVdP(const unsigned int nodeIndex,
        const std::map<std::string,double> & params,
        const std::string & paramName
        ) const
{
  if (paramName.compare("beta") == 0)
    return nodeValues_[nodeIndex];
  else
    return 0.0;
}
// ============================================================================
} // namespace ScalarField
} // namespace Nosh
