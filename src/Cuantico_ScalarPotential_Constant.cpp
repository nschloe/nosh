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

#include "Cuantico_ScalarPotential_Constant.hpp"

#include <Epetra_LocalMap.h>
#include <Epetra_Vector.h>

namespace Cuantico {
namespace ScalarPotential {
// ============================================================================
Constant::
Constant(const double alpha ):
  alpha_( alpha )
{
}
// ============================================================================
Constant::
~Constant()
{
}
// ============================================================================
double
Constant::
getV(const unsigned int nodeIndex,
     const Teuchos::Array<double> & p
     ) const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY(p.length(), 1);
#endif
  return alpha_ + p[0];
}
// ============================================================================
double
Constant::
getdVdP(const unsigned int nodeIndex,
        const unsigned int parameterIndex,
        const Teuchos::Array<double> & p
        ) const
{
  TEUCHOS_ASSERT_EQUALITY(parameterIndex, 0);
  return 1.0;
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<double> >
Constant::
get_p_init() const
{
  return Teuchos::rcp(new Teuchos::Array<double>(1, 0.0));
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
Constant::
get_p_names() const
{
  return Teuchos::rcp(new Teuchos::Array<std::string>(1, "T"));
}
// ============================================================================
} // namespace ScalarPotential
} // namespace Cuantico
