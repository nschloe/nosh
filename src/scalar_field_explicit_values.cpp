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

#include "scalar_field_explicit_values.hpp"

#include <map>
#include <string>

#include "mesh.hpp"

namespace nosh
{
namespace scalar_field
{
// ============================================================================
explicit_values::
explicit_values(const nosh::mesh & mesh,
               const std::string & field_name
             ) :
  node_values_(mesh.get_vector(field_name))
{
}
// ============================================================================
explicit_values::
~explicit_values()
{
}
// ============================================================================
const std::map<std::string, double>
explicit_values::
get_parameters() const
{
  std::map<std::string, double> m;
  m["beta"] = 1.0;
  return m;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
explicit_values::
get_v(const std::map<std::string, double> & params) const
{
  Tpetra::Vector<double,int,int> vals(*node_values_);
  // Scale by "beta"
  vals.scale(params.at("beta"));
  return vals;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
explicit_values::
get_dvdp(const std::map<std::string, double> & params,
        const std::string & param_name
      ) const
{
  (void) params;
  if (param_name.compare("beta") == 0) {
    return *node_values_;
  } else {
    return Tpetra::Vector<double,int,int>(
        node_values_->getMap(),
        true // zero out
        );
  }
}
// ============================================================================
} // namespace scalar_field
} // namespace nosh
