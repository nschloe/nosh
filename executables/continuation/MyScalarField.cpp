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
const std::map<std::string,double>
MyScalarField::
getParameters() const
{
  std::map<std::string,double> m;
  m["tau"] = 0.0;
  return m;
}
// ============================================================================
double
MyScalarField::
getV(const unsigned int nodeIndex,
     const std::map<std::string,double> & params
     ) const
{
  // Get nodal coordinates for nodeIndex.
  std::vector<stk::mesh::Entity*> ownedNodes =
    mesh_->getOwnedNodes();
  const DoubleVector X =
    mesh_->getVectorFieldNonconst(ownedNodes[nodeIndex],
                                  "coordinates", 3);

  // Pick out p["tau"].
  std::map<std::string, double>::const_iterator it = params.find("tau");
  TEUCHOS_ASSERT(it != params.end());
  const double & tau = it->second;

  return -1.0 + tau * (-X[0]*X[0] + X[1]*X[1]);
}
// ============================================================================
double
MyScalarField::
getdVdP(const unsigned int nodeIndex,
        const std::map<std::string,double> & params,
        const std::string & paramName
        ) const
{
  if (paramName.compare("tau") == 0)
  {
    std::vector<stk::mesh::Entity*> ownedNodes =
      mesh_->getOwnedNodes();
    const DoubleVector X =
      mesh_->getVectorFieldNonconst(ownedNodes[nodeIndex],
                                    "coordinates", 3);
    return -X[0]*X[0] + X[1]*X[1];
  }
  else
    return 0.0;
}
// ============================================================================
