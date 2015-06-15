// @HEADER
//
//    Mesh class with compatibility to stk_mesh.
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

//#include <map>
//#include <string>
//#include <algorithm>
//#include <vector>
//
//#include <Teuchos_RCP.hpp>
//
#include <stk_mesh/base/MetaData.hpp>
//#include <stk_mesh/base/BulkData.hpp>
//#include <stk_mesh/base/Entity.hpp>
//// #include <stk_mesh/base/Field.hpp>
//#include <stk_mesh/base/FieldBase.hpp>
//// #include <stk_mesh/base/Comm.hpp> // for comm_mesh_counts
//#include <stk_mesh/base/CreateAdjacentEntities.hpp>
#include <stk_mesh/base/GetEntities.hpp>
//// For parallel_sum:
//#include <stk_mesh/base/FieldParallel.hpp>
//#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
//#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
//#include <stk_mesh/base/CreateFaces.hpp>
//#include <stk_mesh/base/SkinMesh.hpp>
//
//#include <stk_io/IossBridge.hpp>
//#include <Ioss_SubSystem.h>
////#include <stk_io/MeshReadWriteUtils.hpp>
//#include <Ionit_Initializer.h>
//#include <Ioss_IOFactory.h>
//#include <Ioss_Region.h>

//#ifdef HAVE_MPI
//// Rebalance
//#include <stk_rebalance/Rebalance.hpp>
//#include <stk_rebalance_utils/RebalanceUtils.hpp>
//#include <stk_rebalance/Partition.hpp>
//#include <stk_rebalance/ZoltanPartition.hpp>
//#endif

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include "MeshTri.hpp"
#include "MeshTetra.hpp"

namespace Nosh
{
// =============================================================================
std::shared_ptr<Nosh::Mesh>
read(
    const std::string &fileName,
    const int index
    )
{
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

  // ---------------------------------------------------------------------------
  // initialize database communication
  //Ioss::Init::Initializer io;

  // If the file is serial, read it with process 0 and embed it
  // in the multiproc context. Load balancing is done later anyways.
//  MPI_Comm readerComm = //myBulkData->parallel();
//#ifdef HAVE_MPI
//  const bool fileIsSerial = fileName.substr(fileName.find_last_of(".") + 1) == "e";
//  if (fileIsSerial && comm.NumProc() > 1) {
//    // reader process
//    int readerProc[1] = {0};
//
//    // Get the group under mcomm.
//    MPI_Group groupWorld;
//    MPI_Comm_group(mcomm, &groupWorld);
//    // Create the new group.
//    MPI_Group peZero;
//    MPI_Group_incl(groupWorld, 1, readerProc, &peZero);
//    // Create the new communicator.
//    MPI_Comm_create(mcomm, peZero, &readerComm);
//  }
//#endif

  const std::string meshType = "exodusII";
  // By default, the Exodus field separator is '_' such that only fields 'A_*'
  // would be recognized as such.  However, VTK is inconsistent in its putting
  // underscores in field names, such that vectors with three components are
  // typically stored as 'AX', 'AY', 'AZ'.
  // This can be worked around where the Exodus file is created or here.
  // The problem about messing around with the input database here is that the
  // output routines will always "correctly" append an underscore. This then
  // creates a difference between VTK- generated and Trilinos-generated files.
  // Setting removing the underscore from both input and output files is
  // possible in a figure version of Trilinos, c.f.  Greg's mail (2012-09-20):
  // '''
  // Basically, if you do the following prior to creating, it should give you
  // a consistent field separator:
  //
  // meshData->m_property_manager.add(Ioss::Property("FIELD_SUFFIX_SEPARATOR", ""));
  // '''
  //
  // Anyways. Here's some code that removes the underscore from the input
  // database. Keep it commented out for now though.
  // <WORKAROUND CODE START>
  //Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(meshType,
  //                                                fileName,
  //                                                Ioss::READ_MODEL,
  //                                                readerComm);
  //TEUCHOS_TEST_FOR_EXCEPT_MSG(dbi == NULL || !dbi->ok(),
  //                     "ERROR: Could not open database '" << fileName
  //                     << "' of type '" << meshType << "'.");
  //// set the vector field label separator
  //dbi->set_field_separator(0);
  //// create region to feed into meshData
  //Ioss::Region *in_region = new Ioss::Region(dbi, "input_model");
  //meshData->m_input_region = in_region;
  // <WORKAROUND CODE END>

  // This checks the existence of the file, checks to see if we can open it,
  // builds a handle to the region and puts it in mesh_data (in_region), and
  // reads the metaData into metaData.
  ioBroker->add_mesh_database(
      fileName,
      meshType,
      stk::io::READ_MESH
      );
  ioBroker->create_input_mesh();
  //stk::mesh::MetaData& ioBroker->meta_data() = ioBroker->meta_data();
  // ---------------------------------------------------------------------------
  // As of now (2012-03-21) there is no way to determine which fields are
  // actually present in the file. The only thing we can do is to declare them,
  // and check if they're full of zeros in the end.
  //VectorFieldType &coordinatesField =
  //  ioBroker->meta_data().declare_field<VectorFieldType>(
  //      stk::topology::NODE_RANK,
  //      "coordinates"
  //      );
  //stk::mesh::put_field(coordinatesField,
  //                     ioBroker->meta_data().universal_part(),
  //                     numDim);
  //stk::io::set_field_role(coordinatesField,
  //                        Ioss::Field::MESH);

  //IntScalarFieldType &procRankField =
  //  ioBroker->meta_data().declare_field<IntScalarFieldType>(
  //      stk::topology::ELEMENT_RANK,
  //      "proc_rank"
  //      );
  //stk::mesh::put_field(procRankField,
  //                     ioBroker->meta_data().universal_part()
  //                     );
  //stk::io::set_field_role(procRankField,
  //                        Ioss::Field::MESH);

  //// real part
  //ScalarFieldType &psir_field =
  //  ioBroker->meta_data().declare_field<ScalarFieldType>(
  //      stk::topology::NODE_RANK,
  //      "psi_R"
  //      );
  //stk::mesh::put_field(psir_field,
  //                     ioBroker->meta_data().universal_part());
  //stk::io::set_field_role(psir_field,
  //                        Ioss::Field::TRANSIENT);

  //// imaginary part
  //ScalarFieldType &psii_field =
  //  ioBroker->meta_data().declare_field<ScalarFieldType>(
  //      stk::topology::NODE_RANK,
  //      "psi_Z"
  //      );
  //stk::mesh::put_field(psii_field,
  //                     ioBroker->meta_data().universal_part());
  //stk::io::set_field_role(psii_field,
  //                        Ioss::Field::TRANSIENT);

  // Magnetic vector potential.
  // Unconditionally assume that the field is 3D (A_X, A_Y, A_Z) even if the
  // domain is two-dimensional. Eventually, we only need the the projection of
  // the field onto the edges of the mesh, this this might be a bit overkill.
  // Until we can explicitly associate fields with edges, though, keep it this
  // way.
  // Also, declare "A" as Ioss::Field::ATTRIBUTE. This makes sure that the data
  // is written out, but only once (hence "attribute") and not once per step.
  // (This is with trilinos-dev as of July 2012.)
  //VectorFieldType &mvpField =
  //  ioBroker->meta_data().declare_field<VectorFieldType>(
  //      stk::topology::NODE_RANK,
  //      "Arrr"
  //      );
  //stk::mesh::put_field(mvpField,
  //                     ioBroker->meta_data().universal_part(),
  //                     numDim);
  //stk::io::set_field_role(mvpField,
  //                        Ioss::Field::ATTRIBUTE);

  // Note:
  // Thickness and V are actually scalar fields. However, they are both
  // declared VectorFields with 1 component here since otherwise
  // SEACAS throws
  //
  //   p=0: *** Caught standard std::exception of type 'std::runtime_error' :
  //
  //   Expr '!(restr->not_equal_stride(tmp))' eval'd to true, throwing.
  //   Error occured at: stk_mesh/stk_mesh/baseImpl/FieldBaseImpl.cpp:240
  //   Error: std::mesh::MetaData::declare_field_restriction FAILED for FieldBaseImpl<double,Cartesian3d>[ name = "thickness" , #states = 1 ] { entity_rank(0) part({UNIVERSAL}) : 1 } WITH INCOMPATIBLE REDECLARATION { entity_rank(0) part({UNIVERSAL}) : 19754528 }
  //
  //// Thickness field. Same as above.
  //VectorFieldType &thicknessField =
  //  ioBroker->meta_data().declare_field<VectorFieldType>(
  //      stk::topology::NODE_RANK,
  //      "thickness"
  //      );
  //stk::mesh::put_field(thicknessField,
  //                     ioBroker->meta_data().universal_part(),
  //                     1);
  //stk::io::set_field_role(thicknessField, Ioss::Field::ATTRIBUTE);

  //// Potential field. Same as above.
  //VectorFieldType &potentialField =
  //  ioBroker->meta_data().declare_field<VectorFieldType>(
  //      stk::topology::NODE_RANK,
  //      "V"
  //      );
  //stk::mesh::put_field(potentialField,
  //                     ioBroker->meta_data().universal_part(),
  //                     1);
  //stk::io::set_field_role(potentialField, Ioss::Field::ATTRIBUTE);
  // ---------------------------------------------------------------------------

  // Read all fields from the input file
  ioBroker->add_all_mesh_fields_as_input_fields();

  // define_input_fields() doesn't like the ATTRIBUTE fields; disable.
  // What was it good for anyways?
//  stk::io::define_input_fields(*meshData,
//                                metaData
//                             );

//  stk::io::put_io_part_attribute(metaData.universal_part());

  // Finalize the setup.
  //ioBroker->meta_data().commit();

  //myBulkData->modification_begin();
//#ifdef HAVE_MPI
//  if (fileIsSerial && comm.NumProc() > 1) {
//    bulkData->modification_begin();
//    Ioss::Region *region = meshData->m_input_region;
//    if (comm.getRank() == 0) {
//      stk::io::process_mesh_bulk_data(region, *(bulkData));
//      stk::io::input_mesh_fields(region, *(bulkData), index+1);
//    }
//    bulkData->modification_end();
//  } else {
//#endif
    ioBroker->populate_bulk_data();

    // Remember: Indices in STK are 1-based. :/
    ioBroker->read_defined_input_fields(index+1);

    // This should be propagated into stk::io
    //Ioss::Region *region = meshData->m_input_region;
    //const Ioss::ElementBlockContainer& elem_blocks = region->get_element_blocks();
    //// Uncomment to print what fields are in the exodus file
    //Ioss::NameList exo_fld_names;
    //elem_blocks[0]->field_describe(&exo_fld_names);
    //for(std::size_t i = 0; i < exo_fld_names.size(); i++)
    //  std::cout << "Found field \"" << exo_fld_names[i] << "\" in exodus file" << std::endl;

    //myBulkData->modification_end();
//#ifdef HAVE_MPI
//  }

  //// Rebalance.
  //stk::mesh::Selector selector(metaData.universal_part());
  //stk::mesh::Selector owned_selector(metaData.locally_owned_part());
  //double imbalance = stk::rebalance::check_balance(*bulkData,
  //                   NULL,
  //                   metaData.node_rank(),
  //                   &selector);
  //if (imbalance > 1.5) {
  //  //if (comm.getRank() == 0)
  //  //  std::cout << "The imbalance is " << imbalance << ". Rebalance!" << std::endl;
  //  // Zoltan graph-based reblancing.
  //  // http://trilinos.sandia.gov/packages/docs/dev/packages/stk/doc/html/group__stk__rebalance__unit__test__module.html
  //  Teuchos::ParameterList lb_method;
  //  lb_method.set("LOAD BALANCING METHOD", "4");
  //  Teuchos::ParameterList graph;
  //  graph.sublist(stk::rebalance::Zoltan::default_parameters_name()) = lb_method;
  //  stk::rebalance::Zoltan zoltan_partition(mcomm, numDim_, graph);

  //  stk::rebalance::rebalance(*bulkData,
  //                            owned_selector,
  //                            &coordinatesField,
  //                            NULL,
  //                            zoltan_partition);

  //  //imbalance = stk::rebalance::check_balance(bulkData,
  //  //                                          NULL,
  //  //                                          metaData.node_rank(),
  //  //                                          &selector);
  //  //if (comm.getRank() == 0)
  //  //  std::cout << "After rebalancing, the imbalance is " << imbalance << "." << std::endl;
  //}
//#endif

//#ifndef NDEBUG
//  // Assert that all processes own nodes
//  std::vector<stk::mesh::Entity> on = buildOwnedNodes_(ioBroker->bulk_data());
//  TEUCHOS_ASSERT_INEQUALITY(on.size(), >, 0);
//#endif

#if 0
  // create edges
  stk::mesh::create_edges(ioBroker->bulk_data());

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
