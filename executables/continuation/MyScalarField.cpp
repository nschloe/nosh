// @HEADER
//
//    Custom scalar potential.
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

#include "MyScalarField.hpp"

#include "Nosh_StkMesh.hpp"

// ============================================================================
MyScalarField::
MyScalarField(const RCP<const Nosh::StkMesh> & mesh):
  mesh_( mesh )
{
}
// ============================================================================
MyScalarField::
~MyScalarField()
{
}
// ============================================================================
RCP<const Array<std::string> >
MyScalarField::
get_p_names() const
{
  return rcp(new Array<std::string>(1, "tau"));
}
// ============================================================================
RCP<const Array<double> >
MyScalarField::
get_p_init() const
{
  return rcp(new Array<double>(1, 0.0));
}
// ============================================================================
double
MyScalarField::
getV(const unsigned int nodeIndex,
     const Array<double> & p
     ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(p.length(), 1);
#endif
  // Get nodal coordinates for nodeIndex.
  std::vector<stk::mesh::Entity*> ownedNodes =
    mesh_->getOwnedNodes();
  const DoubleVector X =
    mesh_->getVectorFieldNonconst(ownedNodes[nodeIndex],
                                  "coordinates", 3);

  return -1.0 + p[0] * (-X[0]*X[0] + X[1]*X[1]);
}
// ============================================================================
double
MyScalarField::
getdVdP(const unsigned int nodeIndex,
        const unsigned int parameterIndex,
        const Array<double> & p
        ) const
{
  TEUCHOS_ASSERT_EQUALITY(parameterIndex, 0);
  return 1.0;
}
// ============================================================================
