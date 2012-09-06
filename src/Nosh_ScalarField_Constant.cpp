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

#include "Nosh_ScalarField_Constant.hpp"

namespace Nosh {
namespace ScalarField {
// ============================================================================
Constant::
Constant(const double c,
         const std::string & param1Name,
         const double param1InitValue
         ):
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
Teuchos::RCP<const Teuchos::Array<std::string> >
Constant::
get_p_names() const
{
  return Teuchos::rcp(new Teuchos::Array<std::string>(1, param1Name_));
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<double> >
Constant::
get_p_init() const
{
  return Teuchos::rcp(new Teuchos::Array<double>(1, param1InitValue_));
}
// ============================================================================
double
Constant::
getV(const unsigned int nodeIndex,
     const std::map<std::string,double> & params
     ) const
{
  double val = c_;

  if (!param1Name_.empty())
  {
    // If a parameter name was given, add its value to c_;
    std::map<std::string, double>::const_iterator it = params.find(param1Name_);
    TEUCHOS_ASSERT(it != params.end());
    val += it->second;
  }

  return val;
}
// ============================================================================
double
Constant::
getdVdP(const unsigned int nodeIndex,
        const std::map<std::string,double> & params,
        const std::string & paramName
        ) const
{
  if (!param1Name_.empty() && paramName.compare(param1Name_) == 0)
    return 1.0;
  else
    return 0.0;
}
// ============================================================================
} // namespace ScalarField
} // namespace Nosh
