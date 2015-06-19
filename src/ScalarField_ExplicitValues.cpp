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

#include "ScalarField_ExplicitValues.hpp"

#include <map>
#include <string>

#include "Mesh.hpp"

namespace Nosh
{
namespace ScalarField
{
// ============================================================================
ExplicitValues::
ExplicitValues(const Nosh::Mesh & mesh,
               const std::string & fieldName
             ) :
  nodeValues_(mesh.getVector(fieldName))
{
}
// ============================================================================
ExplicitValues::
~ExplicitValues()
{
}
// ============================================================================
const std::map<std::string, double>
ExplicitValues::
getParameters() const
{
  std::map<std::string, double> m;
  m["beta"] = 1.0;
  return m;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
ExplicitValues::
getV(const std::map<std::string, double> & params) const
{
  Tpetra::Vector<double,int,int> vals(*nodeValues_);
  // Scale by "beta"
  vals.scale(params.at("beta"));
  return vals;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
ExplicitValues::
getdVdP(const std::map<std::string, double> & params,
        const std::string & paramName
      ) const
{
  (void) params;
  if (paramName.compare("beta") == 0) {
    return *nodeValues_;
  } else {
    return Tpetra::Vector<double,int,int>(
        nodeValues_->getMap(),
        true // zero out
        );
  }
}
// ============================================================================
} // namespace ScalarField
} // namespace Nosh
