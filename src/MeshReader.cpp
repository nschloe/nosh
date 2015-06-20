// @HEADER
//
//    STK mesh reader.
//    Copyright (C) 2015  Nico Schlömer
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
// =============================================================================
#include "MeshReader.hpp"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
//// For parallel_sum:
//#include <stk_mesh/base/FieldParallel.hpp>
//#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CreateEdges.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include "MeshTri.hpp"
#include "MeshTetra.hpp"

namespace Nosh
{
std::shared_ptr<Nosh::Mesh>
read(
    const std::string & fileName,
    const std::set<std::string> & fields,
    const int index
    )
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const auto fillTime =
    Teuchos::TimeMonitor::getNewTimer("Nosh: read()");
  Teuchos::TimeMonitor tm(*fillTime);
#endif

  auto comm = Teuchos::get_shared_ptr(Teuchos::DefaultComm<int>::getComm());

  auto ioBroker = std::make_shared<stk::io::StkMeshIoBroker>(
#ifdef HAVE_MPI
        *(Teuchos::dyn_cast<const Teuchos::MpiComm<int>>(*comm)
        .getRawMpiComm())
#else
        1
#endif
      );

  // How to split the file for mulitproc I/O
  ioBroker->property_add(
      Ioss::Property("DECOMPOSITION_METHOD", "rcb")
      );
  // Take two different fields with one component
  // instead of one field with two components. This works around
  // Ioss's inability to properly read psi_R, psi_Z as a complex variable.
  // (It can handle data_X, data_Y, data_Z though.)
  //const unsigned int neq = 1;

  const std::string meshType = "exodusII";

  // This checks the existence of the file, checks to see if we can open it,
  // builds a handle to the region and puts it in mesh_data (in_region), and
  // reads the metaData into metaData.
  ioBroker->add_mesh_database(
      fileName,
      meshType,
      stk::io::READ_MESH
      );
  ioBroker->create_input_mesh();

  if (fields.size() > 0) {
    // declare some extra fields so we can write them later
    for (const auto & fieldName: fields) {
      auto & metaData = ioBroker->meta_data();
      ScalarFieldType &field = metaData.declare_field<ScalarFieldType>(
          stk::topology::NODE_RANK,
          fieldName
          );
      stk::mesh::put_field(field, metaData.universal_part());
      stk::io::set_field_role(field, Ioss::Field::TRANSIENT);
    }
  } else {
    // Read all fields from the input file
    ioBroker->add_all_mesh_fields_as_input_fields();
  }

  ioBroker->populate_bulk_data();

  // Remember: Indices in STK are 1-based. :/
  ioBroker->read_defined_input_fields(index+1);

  // create edges
  stk::mesh::create_edges(ioBroker->bulk_data());

#if 0
  std::vector<size_t> mesh_counts;
  stk::mesh::comm_mesh_counts(ioBroker->bulk_data(), mesh_counts);

  std::vector<size_t> counts;
  std::vector<size_t> minCounts;
  std::vector<size_t> maxCounts;
  stk::mesh::comm_mesh_counts(ioBroker->bulk_data(), counts, minCounts, maxCounts);
  if (comm->getRank() == 0) {
    std::cout << "===========================" << std::endl;
    std::cout << "nodes,    " << std::setw(10) << counts[0] << ", "
      << "min/max: " << minCounts[0] << "/" << maxCounts[0] << "\n"
      << "edges,    " << std::setw(10) << counts[1] << ", "
      << "min/max: " << minCounts[1] << "/" << maxCounts[1] << "\n"
      << "faces,    " << std::setw(10) << counts[2] << ", "
      << "min/max: " << minCounts[2] << "/" << maxCounts[2] << "\n"
      << "elements, " << std::setw(10) << counts[3] << ", "
      << "min/max: " << minCounts[3] << "/" << maxCounts[3] << "\n"
      << "==========================="
      << std::endl;
    //if (create_edges) {
    //  >~std::cout<< "num edges/second: " << mesh_counts[stk::topology::EDGE_RANK]/all_edge_time << "\t" << all_edge_time << std::end
    //    l;
    //}
    //if (create_faces) {
    //  >~std::cout<< "num faces/second: " << mesh_counts[stk::topology::FACE_RANK]/all_face_time << "\t" << all_face_time << std::end
    //    l;
    //}
  }

  //// add skin
  //stk::mesh::Part & skin_part =
  //  ioBroker->meta_data().declare_part("skin_part");
  //stk::mesh::PartVector add_parts(1, &skin_part);
  //stk::mesh::skin_mesh(ioBroker->bulk_data(), add_parts);
#endif

  // Check if it's a triangular or tetrahedral mesh.
  // Determine the kind of mesh by the first cell on the first process
  int nodesPerCell;
  if (comm->getRank() == 0) {
    // get owned cells
    stk::mesh::Selector select_owned_in_part =
      stk::mesh::Selector(ioBroker->bulk_data().mesh_meta_data().universal_part())
      & stk::mesh::Selector(ioBroker->bulk_data().mesh_meta_data().locally_owned_part());
    std::vector<stk::mesh::Entity> cells;
    stk::mesh::get_selected_entities(
        select_owned_in_part,
        ioBroker->bulk_data().buckets(stk::topology::ELEMENT_RANK),
        cells
        );
#ifndef NDEBUG
    TEUCHOS_ASSERT_INEQUALITY(cells.size(), >, 0);
#endif
    nodesPerCell = ioBroker->bulk_data().num_nodes(cells[0]);
  }
  Teuchos::broadcast(*comm, 0, 1, &nodesPerCell);

  switch (nodesPerCell) {
    case 3:
      return std::make_shared<Nosh::MeshTri>(comm, ioBroker);
    case 4:
      return std::make_shared<Nosh::MeshTetra>(comm, ioBroker);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          "Illegal cell type."
          );
  }
}

}  // namespace Nosh
