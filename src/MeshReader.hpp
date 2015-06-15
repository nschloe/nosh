// @HEADER
//
//    Mesh class with compatibility to stk_mesh.
//    Copyright (C) 2015 Nico Schlömer
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
#ifndef NOSH_MESHREADER_HPP
#define NOSH_MESHREADER_HPP
// =============================================================================
//// includes
//#include <string>
//#include <vector>
//#include <tuple>
//#include <set>
#include <memory>
//
//#include <Teuchos_RCP.hpp>
//#include <Teuchos_DefaultComm.hpp>
//#include <Teuchos_RCPStdSharedPtrConversions.hpp>
//#ifdef NOSH_TEUCHOS_TIME_MONITOR
//#include <Teuchos_Time.hpp>
//#endif
//#include <Tpetra_Vector.hpp>
//#include <Tpetra_CrsGraph.hpp>
//
//#include <stk_mesh/base/Entity.hpp>
//#include <stk_mesh/base/CoordinateSystems.hpp>
//#include <stk_mesh/base/MetaData.hpp>
//#include <stk_io/StkMeshIoBroker.hpp>
//#include <stk_mesh/base/FieldTraits.hpp>
//
//#include <Eigen/Dense>
//
//// forward declarations
//namespace stk
//{
//namespace mesh
//{
//class BulkData;
//}
//} // namespace stk
//
//// typedefs
//typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
//typedef stk::mesh::Field<double> ScalarFieldType;
//typedef stk::mesh::Field<int> IntScalarFieldType;
//typedef std::tuple<stk::mesh::Entity, stk::mesh::Entity> edge;

#include "Mesh.hpp"

namespace Nosh
{

std::shared_ptr<Nosh::Mesh>
read(
    const std::string & fileName,
    const int index = 0
    );

} // namespace Nosh
// =============================================================================
#endif // NOSH_MESHREADER_HPP
