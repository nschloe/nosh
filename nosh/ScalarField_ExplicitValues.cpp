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

#include "nosh/ScalarField_ExplicitValues.hpp"

#include <map>
#include <string>

#include <Epetra_Vector.h>

#include "nosh/StkMesh.hpp"

namespace Nosh
{
namespace ScalarField
{
// ============================================================================
ExplicitValues::
ExplicitValues(const Nosh::StkMesh & mesh,
               const std::string & fieldName
             ) :
  nodeValues_(mesh.createVector(fieldName))
{
}
// ============================================================================
ExplicitValues::
~ExplicitValues()
{
}
// ============================================================================
const std::map<std::string,double>
ExplicitValues::
getInitialParameters() const
{
  std::map<std::string,double> m;
  m["beta"] = 1.0;
  return m;
}
// ============================================================================
const Epetra_Vector
ExplicitValues::
getV(const std::map<std::string,double> & params) const
{
  Epetra_Vector vals(*nodeValues_);

  // Find the value of "beta" and use it as a factor.
  std::map<std::string, double>::const_iterator it = params.find("beta");
  TEUCHOS_ASSERT(it != params.end());
  vals.Scale(it->second);

  return vals;
}
// ============================================================================
const Epetra_Vector
ExplicitValues::
getdVdP(const std::map<std::string,double> & params,
        const std::string & paramName
      ) const
{
  (void) params;
  if (paramName.compare("beta") == 0)
    return *nodeValues_;
  else
    return Epetra_Vector(nodeValues_->Map(), true); // 0.0 overall
}
// ============================================================================
} // namespace ScalarField
} // namespace Nosh
