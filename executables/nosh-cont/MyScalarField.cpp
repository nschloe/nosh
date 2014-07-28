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

#include <map>
#include <string>
#include <vector>

#include "nosh/StkMesh.hpp"

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

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
const std::map<std::string, double>
MyScalarField::
getInitialParameters() const
{
  std::map<std::string, double> m;
  m["tau"] = 0.0;
  return m;
}
// ============================================================================
const Epetra_Vector
MyScalarField::
getV(const std::map<std::string, double> & params) const
{
  // Pick out p["tau"].
  std::map<std::string, double>::const_iterator it = params.find("tau");
  TEUCHOS_ASSERT(it != params.end());
  const double & tau = it->second;

  std::vector<stk_classic::mesh::Entity*> ownedNodes =
    mesh_->getOwnedNodes();

  Epetra_Vector vals(*(mesh_->getNodesMap()));

  for (unsigned int k = 0; k < ownedNodes.size(); k++) {
    // Get nodal coordinates.
    const DoubleVector X =
      mesh_->getVectorFieldNonconst(ownedNodes[k],
                                    "coordinates", 3);
    vals[k] = -1.0 + tau * (-X[0]*X[0] + X[1]*X[1]);
  }

  return vals;
}
// ============================================================================
const Epetra_Vector
MyScalarField::
getdVdP(const std::map<std::string, double> & params,
        const std::string & paramName
       ) const
{
  // Silence warning about unused params.
  (void) params;

  // Create vals as zeroed-out vector.
  Epetra_Vector vals(*(mesh_->getNodesMap()));

  if (paramName.compare("tau") == 0) {
    std::vector<stk_classic::mesh::Entity*> ownedNodes =
      mesh_->getOwnedNodes();

    for (unsigned int k = 0; k < ownedNodes.size(); k++) {
      // Get nodal coordinates.
      const DoubleVector X =
        mesh_->getVectorFieldNonconst(ownedNodes[k],
                                      "coordinates", 3);
      vals[k] = -X[0]*X[0] + X[1]*X[1];
    }
  }

  return vals;
}
// ============================================================================
