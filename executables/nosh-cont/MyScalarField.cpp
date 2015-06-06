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

#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include "nosh/StkMesh.hpp"

// ============================================================================
MyScalarField::
MyScalarField(const std::shared_ptr<const Nosh::StkMesh> & mesh):
  mesh_(mesh)
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
getParameters() const
{
  std::map<std::string, double> m;
  m["tau"] = 0.0;
  return m;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
MyScalarField::
getV(const std::map<std::string, double> & params) const
{
  const double tau = params.at("tau");

  auto ownedNodes = mesh_->getOwnedNodes();

  Tpetra::Vector<double,int,int> vals(Teuchos::rcp(mesh_->getNodesMap()));
  auto vData = vals.getDataNonConst();

  const VectorFieldType & coordsField = mesh_->getNodeField("coordinates");

  for (size_t k = 0; k < ownedNodes.size(); k++) {
    // Get nodal coordinates.
    const Eigen::Vector3d X = mesh_->getNodeValue(coordsField, ownedNodes[k]);
    vData[k] = -1.0 + tau * (-X[0]*X[0] + X[1]*X[1]);
  }

  return vals;
}
// ============================================================================
const Tpetra::Vector<double,int,int>
MyScalarField::
getdVdP(
    const std::map<std::string, double> & params,
    const std::string & paramName
    ) const
{
  // Silence warning about unused params.
  (void) params;

  // Create vals as zeroed-out vector.
  Tpetra::Vector<double,int,int> vals(Teuchos::rcp(mesh_->getNodesMap()), true);

  if (paramName.compare("tau") == 0) {
    std::vector<stk::mesh::Entity> ownedNodes =
      mesh_->getOwnedNodes();

    const VectorFieldType & coordsField = mesh_->getNodeField("coordinates");

    auto vData = vals.getDataNonConst();

    for (size_t k = 0; k < ownedNodes.size(); k++) {
      // Get nodal coordinates.
      const Eigen::Vector3d X = mesh_->getNodeValue(coordsField, ownedNodes[k]);
      vData[k] = -X[0]*X[0] + X[1]*X[1];
    }
  }

  return vals;
}
// ============================================================================
