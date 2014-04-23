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

#include "nosh/ScalarField_Constant.hpp"

#include <map>
#include <string>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

namespace Nosh
{
namespace ScalarField
{
// ============================================================================
Constant::
Constant(const Nosh::StkMesh & mesh,
         const double c,
         const std::string & param1Name,
         const double param1InitValue
        ):
  map_(mesh.getNodesMap()),
  c_( c ),
  param1Name_ (param1Name),
  param1InitValue_( param1InitValue )
{
}
// ============================================================================
Constant::
~Constant()
{
}
// ============================================================================
const std::map<std::string,double>
Constant::
getInitialParameters() const
{
  std::map<std::string,double> m;
  if (!param1Name_.empty())
    m[param1Name_] = param1InitValue_;
  return m;
}
// ============================================================================
const Epetra_Vector
Constant::
getV(const std::map<std::string,double> & params) const
{
  // Create constant-valued vector.
  Epetra_Vector vals(*map_);

  std::map<std::string, double>::const_iterator it = params.find(param1Name_);
  if (it != params.end())
    vals.PutScalar(c_ + it->second);
  else
    vals.PutScalar(c_);

  return vals;
}
// ============================================================================
const Epetra_Vector
Constant::
getdVdP(const std::map<std::string,double> & params,
        const std::string & paramName
       ) const
{
  (void) params;
  // Create zeroed-out vector.
  Epetra_Vector vals(*map_, true);
  if (paramName.compare(param1Name_) == 0)
    vals.PutScalar(1.0);

  return vals;
}
// ============================================================================
} // namespace ScalarField
} // namespace Nosh
