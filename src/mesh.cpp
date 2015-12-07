// @HEADER
//
//    Mesh class with compatibility to stk_mesh.
//    Copyright (C) 2010--2012  Nico Schlömer
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
#include "mesh.hpp"

#include <moab/Core.hpp>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#  include <Teuchos_TimeMonitor.hpp>
#endif

namespace nosh
{
// =============================================================================
mesh::
mesh(
    const std::shared_ptr<const Teuchos::Comm<int>> & _comm,
    const std::shared_ptr<moab::ParallelComm> & mcomm,
    const std::shared_ptr<moab::Core> & mb
    ) :
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  write_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: mesh::write")),
  multi_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: mesh::get_multi_vector")),
#endif
  comm(_comm),
  mb_(mb),
  mcomm_(mcomm),
  //owned_nodes_(this->build_owned_nodes_()),
  nodes_map_(this->get_map_(this->get_owned_gids_())),
  nodes_overlap_map_(this->get_map_(this->get_overlap_gids_())),
  complex_map_(this->get_map_(this->complexify_(this->get_owned_gids_()))),
  complex_overlap_map_(
    this->get_map_(this->complexify_(this->get_overlap_gids_()))
    )
  ,edge_data_(this->build_edge_data_())
  ,edge_lids(build_edge_lids_())
  ,edge_lids_complex(build_edge_lids_complex_())
{
// TODO resurrect
//#ifndef NDEBUG
//  // Assert that all processes own nodes
//  TEUCHOS_ASSERT_INEQUALITY(owned_nodes_.size(), >, 0);
//#endif
  std::cout << "mesh::mesh" << std::endl;
}
// =============================================================================
mesh::
~mesh()
{
}
// =============================================================================
//std::shared_ptr<Tpetra::Vector<double,int,int>>
//mesh::
//field_to_vector_(const ScalarFieldType &field) const
//{
//  // Get overlap nodes.
//  const auto & overlap_nodes = this->get_overlap_nodes();
//
//  // Create vector with this respective map.
//  auto vector = std::make_shared<Tpetra::Vector<double,int,int>>(
//      Teuchos::rcp(this->overlap_map())
//      );
//
//  auto v_data = vector->getDataNonConst();
//
//#ifndef NDEBUG
//  TEUCHOS_ASSERT_EQUALITY(v_data.size(), 2*overlap_nodes.size());
//#endif
//
//  // Fill the vector with data from the file.
//  for (unsigned int k = 0; k < overlap_nodes.size(); k++) {
//    double* vals = stk::mesh::field_data(field, overlap_nodes[k]);
//    // Check if the field is actually there.
//#ifndef NDEBUG
//    TEUCHOS_ASSERT(vals != NULL);
//    //*out << "WARNING: _value for node " << k << " not found.\n" <<
//    //  "Probably there is no field given with the state. Using default."
//    //  << std::endl;
//#endif
//    v_data[k] = vals[0];
//  }
//
//#ifndef NDEBUG
//  // Use NormInf as it's robust against overlapping maps.
//  const double r = vector->normInf();
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(
//      r != r || r > 1.0e100,
//      "The input data seems flawed. Abort."
//      );
//#endif
//
//  return vector;
//}
// =============================================================================
//std::shared_ptr<Tpetra::MultiVector<double,int,int>>
//mesh::
//field_to_vector_(
//    const vector_fieldType &field,
//    const int num_components
//    ) const
//{
//  // Get overlap nodes.
//  const auto & overlap_nodes = this->get_overlap_nodes();
//
//  // Create vector with this respective map.
//  auto vector = std::make_shared<Tpetra::MultiVector<double,int,int>>(
//      Teuchos::rcp(this->overlap_map()),
//      num_components
//      );
//
//  std::vector<Teuchos::ArrayRCP<double>> data(num_components);
//  for (int i = 0; i < num_components; i++) {
//    data[i] = vector->getDataNonConst(i);
//#ifndef NDEBUG
//    TEUCHOS_ASSERT_EQUALITY(data[i].size(), overlap_nodes.size());
//#endif
//  }
//
//  // Fill the vector with data from the file.
//  for (unsigned int k = 0; k < overlap_nodes.size(); k++) {
//    const double * const vals =
//      stk::mesh::field_data(field, overlap_nodes[k]);
//#ifndef NDEBUG
//    // Check if the field is actually there.
//    TEUCHOS_TEST_FOR_EXCEPT_MSG(
//        vals == NULL,
//        "Field value for node " << k << " not found.\n" <<
//        "Probably there is no field given with the state."
//        );
//#endif
//    // Copy over.
//    // A multivector isn't actually a good data structure for this.  What would
//    // be needed is a vector where each entry has k components. This way, the
//    // data could stick together.
//    for (int i = 0; i < num_components; i++) {
//      data[i][k] = vals[i];
//    }
//  }
//
//#ifndef NDEBUG
//  // Check for NaNs and uninitialized data.
//  std::vector<double> r(num_components);
//  // Use NormInf as it's robust against overlapping maps.
//  vector->normInf(Teuchos::ArrayView<double>(r));
//  for (int i = 0; i < num_components; i++) {
//    TEUCHOS_TEST_FOR_EXCEPT_MSG(
//        r[i] != r[i] || r[i] > 1.0e100,
//        "The input data seems flawed. Abort."
//        );
//  }
//#endif
//
//  return vector;
//}
// =============================================================================
void
mesh::
open_file(
    const std::string &output_file
    )
{
#if 0
  output_channel_ = io_broker_->create_output_mesh(
      output_file,
      stk::io::WRITE_RESULTS
      );
  const stk::mesh::FieldVector &fields = io_broker_->meta_data().get_fields();
  for (size_t i=0; i < fields.size(); i++) {
    if (*stk::io::get_field_role(*fields[i]) == Ioss::Field::TRANSIENT) {
      io_broker_->add_field(output_channel_, *fields[i]);
    }
  }
#endif
  return;
}
// =============================================================================
void
mesh::
write(const std::string & filename) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*write_time_);
#endif

  this->mb_->write_mesh(filename.c_str());

  return;
}
// =============================================================================
std::vector<double>
mesh::
get_data(
  const std::string & tag_name,
  const moab::Range & range
  ) const
{
  moab::ErrorCode ierr;

  moab::Tag tag;
  ierr = mb_->tag_get_handle(tag_name.c_str(), tag);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  int length;
  ierr = mb_->tag_get_length(tag, length);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  moab::DataType type;
  ierr = mb_->tag_get_data_type(tag, type);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  TEUCHOS_ASSERT_EQUALITY(type, moab::DataType::MB_TYPE_DOUBLE);

  const int num_data = length * range.size();
  std::vector<double> data(num_data);
  ierr = mb_->tag_get_data(tag, range, &data[0]);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  return data;
}
// =============================================================================
std::vector<double>
mesh::
get_coords(
  const moab::EntityHandle vertex
  ) const
{
  moab::ErrorCode ierr;

  std::vector<double> coords(3);
  ierr = this->mb_->get_coords(&vertex, 1, &coords[0]);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  return coords;
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh::
get_vector(const std::string & tag_name) const
{
#if 0
  const ScalarFieldType * const field =
    io_broker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        field_name
        );

  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      field == NULL,
      "Scalar field \"" << field_name << "\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  return this->field_to_vector_(*field);
#endif

  return std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(this->overlap_map())
      );
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh::
get_complex_vector(const std::string & tag_name) const
{
  // get data for all vertices
  moab::ErrorCode ierr;
  moab::Range verts;
  ierr = this->mb_->get_entities_by_dimension(0, 0, verts);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  auto data = this->get_data(tag_name, verts);

  TEUCHOS_ASSERT_EQUALITY(
    data.size(),
    this->complex_overlap_map_->getNodeNumElements()
    );

  // Set vector values from an existing array (copy)
  return std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(this->complex_overlap_map_),
      Teuchos::ArrayView<double>(data)
      );
}
// =============================================================================
std::shared_ptr<Tpetra::MultiVector<double,int,int>>
mesh::
get_multi_vector(const std::string & tag_name) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*multi_time_);
#endif
  // get data for all vertices
  moab::ErrorCode ierr;
  moab::Range verts;
  ierr = mb_->get_entities_by_dimension(0, 0, verts);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  auto data = this->get_data(tag_name, verts);

  // MOAB's ordering is
  //   x0, y0, z0, x1, y1, z1, ...
  // However,Tpetra::MultiVector's constructor needs the data ordered like
  //   x0, x1, ..., xn, y0, y1, ...
  // Hence, reorder.
  std::vector<double> new_data(data.size());
  const size_t length = data.size() / verts.size();
  for (size_t i = 0; i < length; i++) {
    for (size_t j = 0; j < verts.size(); j++) {
      new_data[j + i*verts.size()] = data[i + j*length];
    }
  }

  TEUCHOS_ASSERT_EQUALITY(
    data.size(),
    length * this->overlap_map()->getNodeNumElements()
    );

  // Set vector values from an existing array (copy)
  return std::make_shared<Tpetra::MultiVector<double,int,int>>(
      Teuchos::rcp(this->overlap_map()),
      Teuchos::ArrayView<double>(new_data),
      verts.size(),
      length
      );
}
// =============================================================================
void
mesh::
insert_vector(
    const Tpetra::Vector<double,int,int> & x,
    const std::string & field_name
    ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*write_time_);
#endif

#if 0
  ScalarFieldType * x_field =
    io_broker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        field_name
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      x_field == NULL,
      "Scalar field \"" << field_name << "\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  // Zero out all nodal values, including the overlaps.
  const auto & overlap_nodes = this->get_overlap_nodes();
  for (unsigned int k = 0; k < overlap_nodes.size(); k++) {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data(*x_field, overlap_nodes[k]);
    *localPsiR = 0.0;
  }

  auto x_data = x.getData();

  // Set owned nodes.
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(x_data.size(), owned_nodes_.size());
#endif
  for (unsigned int k = 0; k < owned_nodes_.size(); k++) {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data(*x_field, owned_nodes_[k]);
    *localPsiR = x_data[k];
  }

  // This communication updates the field values on un-owned nodes it is
  // correct because the zeroSolutionField above zeros them all and the
  // getSolutionField only sets the owned nodes.
  // TODO combine these fields into a vector of fields
  std::vector<stk::mesh::FieldBase*> tmp(1, x_field);
  stk::mesh::parallel_sum(io_broker_->bulk_data(), tmp);

#endif
  return;
}
// =============================================================================
std::vector<moab::EntityHandle>
mesh::
get_owned_cells() const
{
  std::vector<moab::EntityHandle> cells;
  throw "mesh::get_owned_cells";
#if 0
  // get owned elements
  stk::mesh::Selector select_owned_in_part =
    stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().universal_part())
    & stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().locally_owned_part());
  stk::mesh::get_selected_entities(
      select_owned_in_part,
      io_broker_->bulk_data().buckets(stk::topology::ELEMENT_RANK),
      cells
      );
#endif
  return cells;
}
// =============================================================================
std::vector<moab::EntityHandle>
mesh::
get_overlap_edges() const
{
  std::vector<moab::EntityHandle> edges;
#if 0
  // get overlap edges
  stk::mesh::Selector select_overlap_in_part =
    stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().universal_part())
    & (stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().locally_owned_part())
       |stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().globally_shared_part()));

  stk::mesh::get_selected_entities(
      select_overlap_in_part,
      io_broker_->bulk_data().buckets(stk::topology::EDGE_RANK),
      edges
      );
#endif
  return edges;
}
// =============================================================================
//const vector_fieldType &
//mesh::
//get_node_field(const std::string & field_name) const {
//  const vector_fieldType * const field =
//    io_broker_->bulk_data().mesh_meta_data().get_field<vector_fieldType>(
//        stk::topology::NODE_RANK,
//        field_name
//        );
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(
//      field == NULL,
//      "Invalid field name \"" << field_name << "\"."
//      );
//  return *field;
//}
// =============================================================================
const std::vector<Teuchos::Tuple<int,2>>
mesh::
build_edge_lids_() const
{
  std::cout << ">> build_edge_lids_" << std::endl;
  const std::vector<edge> edges = this->my_edges();

  std::vector<Teuchos::Tuple<int,2>> _edge_lids(edges.size());

  for (std::size_t k = 0; k < edges.size(); k++) {
    _edge_lids[k] = Teuchos::tuple(
        (int)this->local_index(std::get<0>(edges[k])),
        (int)this->local_index(std::get<1>(edges[k]))
        );
  }

  std::cout << "   build_edge_lids_ >>" << std::endl;
  return _edge_lids;
}
// =============================================================================
const std::vector<Teuchos::Tuple<int,4>>
mesh::
build_edge_lids_complex_() const
{
  const std::vector<edge> edges = this->my_edges();

  std::vector<Teuchos::Tuple<int,4>> _edge_lids_complex(edges.size());

  for (std::size_t k = 0; k < edges.size(); k++) {
    int lidT0 = this->local_index(std::get<0>(edges[k]));
    int lidT1 = this->local_index(std::get<1>(edges[k]));
    _edge_lids_complex[k] =
      Teuchos::tuple(
          2*lidT0, 2*lidT0 + 1,
          2*lidT1, 2*lidT1 + 1
          );
  }

  return _edge_lids_complex;
}
// =============================================================================
const std::vector<int>
mesh::
get_owned_gids_() const
{
  moab::ErrorCode ierr;
  const auto mb = this->mcomm_->get_moab();

  // get owned entities
  moab::Range verts;
  ierr = mb->get_entities_by_dimension(0, 0, verts);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  // get the corresponding global IDs
  moab::Tag gid;
  ierr = mb->tag_get_handle("GLOBAL_ID", gid);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  std::vector<int> global_ids(verts.size());
  ierr = mb->tag_get_data(gid, verts, &global_ids[0]);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  return global_ids;
}
// =============================================================================
const std::vector<int>
mesh::
get_overlap_gids_() const
{
  moab::ErrorCode ierr;
  const auto mb = this->mcomm_->get_moab();

  // get owned
  moab::Range all;
  ierr = mb->get_entities_by_dimension(0, 0, all);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  // Get entities shared with all other processors
  moab::Range shared;
  ierr = mcomm_->get_shared_entities(-1, shared, 0);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  // merge
  all.merge(shared);

  // get the corresponding global IDs
  moab::Tag gid;
  ierr = mb->tag_get_handle("GLOBAL_ID", gid);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  std::vector<int> global_ids(all.size());
  ierr = mb->tag_get_data(gid, all, &global_ids[0]);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  // std::cout << "rank " << comm->getRank() << " " << global_ids.size() << std::endl;
  // for (size_t k = 0; k < global_ids.size(); k++) {
  //   std::cout << "  " << comm->getRank() << " " << global_ids[k] << std::endl;
  // }

  return global_ids;
}
// =============================================================================
const std::vector<int>
mesh::
complexify_(const std::vector<int> & ids) const
{
  std::vector<int> complex_ids(2 * ids.size());
  for (size_t k=0; k < ids.size(); k++) {
    complex_ids[2*k] = 2 * ids[k];
    complex_ids[2*k+1] = 2 * ids[k] + 1;
  }

  return complex_ids;
}
// =============================================================================
std::shared_ptr<Tpetra::Map<int,int>>
mesh::
get_map_(const std::vector<int> & ids) const
{
  return std::make_shared<Tpetra::Map<int,int>>(
      -1,
      ids,
      0,
      Teuchos::rcp(this->comm)
      );
}
// =============================================================================
mesh::edges_container
mesh::
build_edge_data_()
{
  moab::ErrorCode rval;

  // get the number of 3D entities
  int num3d = 0;
  rval = this->mb_->get_number_entities_by_dimension(0, 3, num3d);
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);

  const int dim = (num3d > 0) ? 3 : 2;

  std::cout << "dim = " << dim << std::endl;

  // Get regions, by dimension, so we stay generic to entity type
  moab::Range elems;
  rval = this->mb_->get_entities_by_dimension(0, dim, elems);
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
  std::cout << "Number of elements is " << elems.size() << std::endl;
  //for (size_t k = 0; k < elems.size(); k++) {
  //  std::cout << "elems[" << k << "] = " << elems[k] << std::endl;
  //}

  // get and create all edges adjacent to cells
  moab::Range edges;
  rval = this->mb_->get_adjacencies(
      elems,
      1,
      true,
      edges,
      moab::Interface::UNION
      );
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
  std::cout << "Number of edges adjacent to cells: " << edges.size() << std::endl;

  // for testing: TODO remove
  moab::Range verts2;
  rval = mb_->get_entities_by_dimension(0, 0, verts2);
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
  for (size_t k = 0; k < verts2.size(); k++) {
    auto tmp = verts2[k];
    std::vector<double> coords(3);
    rval = this->mb_->get_coords(&tmp, 1, &coords[0]);
    TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
    std::cout << "vert " << k << "   coords"
      << " " << coords[0]
      << " " << coords[1]
      << " " << coords[2]
      << std::endl;
  }
  // std::vector<double> coords;
  // ierr = mb_->get_vertex_coordinates(coords);
  // TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  // create cell->edge relation
  std::vector<std::vector<moab::EntityHandle>> cell_edges(elems.size());
  for (size_t k = 0; k < elems.size(); k++) {
    // TODO don't use tmp
    moab::EntityHandle tmp = elems[k];
    rval = this->mb_->get_adjacencies(
        &tmp, 1,
        1,
        true,
        cell_edges[k],
        moab::Interface::UNION
        );
    TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);

    //std::cout << "cell_edges[" << k << "] ="
    //  << " " << cell_edges[k][0]
    //  << " " << cell_edges[k][1]
    //  << " " << cell_edges[k][2]
    //  // << " " << out[0]
    //  // << " " << out[1]
    //  // << " " << out[2]
    //  << std::endl;
  }

  // create edge->node relation
  std::vector<std::tuple<moab::EntityHandle, moab::EntityHandle>>
    edge_nodes(edges.size());
  for (size_t k = 0; k < edges.size(); k++) {
    moab::Range verts;
    // TODO don't use tmp
    std::vector<moab::EntityHandle> tmp(1);
    tmp[0] = edges[k];
    rval = this->mb_->get_adjacencies(
        &tmp[0], 1,
        0,
        true,
        verts,
        moab::Interface::UNION
        );
    TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
    edge_nodes[k] = std::make_tuple(verts[0], verts[1]);
    std::cout << "edge " << k << ", verts " << verts[0] << " " << verts[1] << std::endl;
  }

  mesh::edges_container edge_data = {edge_nodes, cell_edges};

  // for testing
  moab::Range verts;
  rval = this->mb_->get_adjacencies(
      elems,
      0,
      true,
      verts,
      moab::Interface::UNION
      );
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);
  std::cout << "Number of verts adjacent to cells: " << verts.size() << std::endl;
  // for (size_t k = 0; k < verts.size(); k++) {
  //   std::cout << "vert " << k << " " << verts[k] << std::endl;
  // }

#if 0
  std::vector<moab::EntityHandle> cells = this->get_owned_cells();
  size_t num_local_cells = cells.size();

  mesh::edges_container edge_data = {
    // Local edge ID -> Local node IDs.
    std::vector<std::tuple<moab::EntityHandle, moab::EntityHandle> >(),
    // Local cell ID -> Local edge IDs.
    std::vector<std::vector<int>>(num_local_cells)
    };

  // Used to determine if an edge sits on the domain boundary
  std::vector<size_t> num_cells_per_edge();

  // This std::map keeps track of how nodes and edges are connected. If
  // `node_edges((3,4)) == 17`, then the nodes (3,4) are connected by edge 17.
  // Unfortunately, std::tuples can't be compared with '<'. Provide a function
  // pointer that implements lexicographic comparison.
  // See http://www.cplusplus.com/reference/stl/map/map/.
  std::map<std::tuple<moab::EntityHandle, moab::EntityHandle>, int> nodesEdge;

  //const EntityComp ec(io_broker_->bulk_data());

  // Loop over all owned cells.
  unsigned int edge_local_id = 0;
  for (size_t cellLID = 0; cellLID < num_local_cells; cellLID++) {
    // Loop over all pairs of local nodes.
    moab::EntityHandle const * local_nodes
      = io_broker_->bulk_data().begin_nodes(cells[cellLID]);
    const size_t num_local_nodes =
      io_broker_->bulk_data().num_nodes(cells[cellLID]);

    //stk::mesh::PairIterRelation nodesIterator =
    //  cells[cellLID]->relations(meta_data.node_rank());
    //unsigned int num_local_nodes = nodesIterator.size();
    size_t num_local_edges = num_local_nodes * (num_local_nodes - 1) / 2;

    edge_data.cell_edges[cellLID] = std::vector<int>(num_local_edges);

    // Gather the node entities.
    std::vector<moab::EntityHandle> nodes(num_local_nodes);
    for (size_t k = 0; k < num_local_nodes; k++) {
      nodes[k] = local_nodes[k];
    }

    // Sort nodes. This is necessary to make sure that the tuples formed below
    // are always sorted such they are unique keys (and {3,7}, {7,3} are
    // recognized as the same edge).
    std::sort(nodes.begin(), nodes.end());

    // In a simplex, the edges are exactly the connection between each pair of
    // nodes. Hence, loop over pairs of nodes.
    unsigned int edge_index = 0;
    edge edge_nodes;
    for (size_t e0 = 0; e0 < num_local_nodes; e0++) {
      std::get<0>(edge_nodes) = nodes[e0];
      for (size_t e1 = e0+1; e1 < num_local_nodes; e1++) {
        std::get<1>(edge_nodes) = nodes[e1];
        // As nodes are sorted and by their identifiers, edge_nodes are sorted
        // too. This is necessary as otherwise the edge {3,7} could not be
        // identified as {7,3}.
        // Check if edge_nodes is in the map.
        auto it = nodesEdge.find(edge_nodes);
        if (it != nodesEdge.end()) {
          // Edge is already accounted for.
          edge_data.cell_edges[cellLID][edge_index] = it->second;
        } else {  // Edge not found -- insert it.
          nodesEdge[edge_nodes] = edge_local_id; // for householding in this method
          edge_data.edge_nodes.push_back(edge_nodes); // for looping over edges
          edge_data.cell_edges[cellLID][edge_index] = edge_local_id; // for this->compute_edge_coefficients_
          edge_local_id++;
        }
        edge_index++;
      }
    }
  }
#endif
  return edge_data;
}
// =============================================================================
Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
mesh::
build_graph() const
{
  std::cout << ">> build_graph" << std::endl;
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from
  // the owned_nodes map, and the overlapping one from the
  // overlap_nodes.
  // Let's illustrate the implications with the example of the matrix
  //   [ 2 1   ]
  //   [ 1 2 1 ]
  //   [   1 2 ].
  // Suppose subdomain 1 consists of node 1, subdomain 2 of node 3,
  // and node 2 forms the boundary between them.
  // For two processes, if process 1 owns nodes 1 and 2, the matrix
  // will be split as
  //   [ 2 1   ]   [       ]
  //   [ 1 2 1 ] + [       ]
  //   [       ]   [   1 2 ].
  // The vectors always need to have a unique map (otherwise, norms
  // cannot be computed by Epetra), so let's assume they have the
  // map ([1, 2], [3]).
  // The communucation for a matrix-vector multiplication Ax=y
  // needs to be:
  //
  //   1. Communicate x(3) to process 1.
  //   2. Communicate x(2) to process 2.
  //   3. Compute.
  //
  // If the matrix is split up like
  //   [ 2 1   ]   [       ]
  //   [ 1 1   ] + [   1 1 ]
  //   [       ]   [   1 2 ]
  // (like the overlap map suggests), then any Ax=y comes down to:
  //
  //   1. Communicate x(2) to process 2.
  //   2. Compute.
  //   3. Communicate (part of) y(2) to process 1.
  //
  // In the general case, assuming that the number of nodes adjacent
  // to a boundary (on one side) are approximately the number of
  // nodes on that boundary, there is not much difference in
  // communication between the patterns.
  // What does differ, though, is the workload on the processes
  // during the computation phase: Process 1 that owns the whole
  // boundary, has to compute more than process 2.
  // Notice, however, that the total number of computations is
  // lower in scenario 1 (7 vs. 8 FLOPs); the same is true for
  // storage.
  // Hence, it comes down to the question whether or not the
  // mesh generator provided a fair share of the boundary nodes.
  // If yes, then scenario 1 will yield approximately even
  // computation times; if not, then scenario 2 will guarantee
  // equal computation times at the cost of higher total
  // storage and computation needs.
  //
  // Remark:
  // This matrix will later be fed into ML. ML has certain restrictions as to
  // what maps can be used. One of those is that RowMatrixRowMap() and
  // getRangeMap must be the same, and, if the matrix is square,
  // getRangeMap and getDomainMap must coincide too.
  //
  const auto nonoverlap_map = this->map();
#ifndef NDEBUG
  TEUCHOS_ASSERT(nonoverlap_map);
#endif

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  auto graph = Teuchos::rcp(new Tpetra::CrsGraph<int, int>(
      Teuchos::rcp(nonoverlap_map),
      Teuchos::rcp(nonoverlap_map),
      0
      ));

  const std::vector<edge> edges = this->my_edges();

  std::cout << "hasColMap " << graph->hasColMap() << std::endl;
  std::cout << "isLocallyIndexed " << graph->isLocallyIndexed() << std::endl;
  std::cout << "isGloballyIndexed " << graph->isGloballyIndexed() << std::endl;

  // Loop over all edges and put entries wherever two nodes are connected.
  // TODO check if we can use LIDs here
  for (size_t k = 0; k < edges.size(); k++) {
    const Teuchos::Tuple<int,2> & idx = this->edge_lids[k];
    for (int i = 0; i < 2; i++) {
      graph->insertLocalIndices(idx[i], idx);
    }
  }
  graph->fillComplete();

  std::cout << "   build_graph >>" << std::endl;
  return graph;
}
// =============================================================================
Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
mesh::
build_complex_graph() const
{
  std::cout << ">> build_complex_graph" << std::endl;
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from
  // the owned_nodes map, and the overlapping one from the
  // overlap_nodes.
  // Let's illustrate the implications with the example of the matrix
  //   [ 2 1   ]
  //   [ 1 2 1 ]
  //   [   1 2 ].
  // Suppose subdomain 1 consists of node 1, subdomain 2 of node 3,
  // and node 2 forms the boundary between them.
  // For two processes, if process 1 owns nodes 1 and 2, the matrix
  // will be split as
  //   [ 2 1   ]   [       ]
  //   [ 1 2 1 ] + [       ]
  //   [       ]   [   1 2 ].
  // The vectors always need to have a unique map (otherwise, norms
  // cannot be computed by Epetra), so let's assume they have the
  // map ([1, 2], [3]).
  // The communucation for a matrix-vector multiplication Ax=y
  // needs to be:
  //
  //   1. Communicate x(3) to process 1.
  //   2. Communicate x(2) to process 2.
  //   3. Compute.
  //
  // If the matrix is split up like
  //   [ 2 1   ]   [       ]
  //   [ 1 1   ] + [   1 1 ]
  //   [       ]   [   1 2 ]
  // (like the overlap map suggests), then any Ax=y comes down to:
  //
  //   1. Communicate x(2) to process 2.
  //   2. Compute.
  //   3. Communicate (part of) y(2) to process 1.
  //
  // In the general case, assuming that the number of nodes adjacent
  // to a boundary (on one side) are approximately the number of
  // nodes on that boundary, there is not much difference in
  // communication between the patterns.
  // What does differ, though, is the workload on the processes
  // during the computation phase: Process 1 that owns the whole
  // boundary, has to compute more than process 2.
  // Notice, however, that the total number of computations is
  // lower in scenario 1 (7 vs. 8 FLOPs); the same is true for
  // storage.
  // Hence, it comes down to the question whether or not the
  // mesh generator provided a fair share of the boundary nodes.
  // If yes, then scenario 1 will yield approximately even
  // computation times; if not, then scenario 2 will guarantee
  // equal computation times at the cost of higher total
  // storage and computation needs.
  //
  // Remark:
  // This matrix will later be fed into ML. ML has certain restrictions as to
  // what maps can be used. One of those is that RowMatrixRowMap() and
  // getRangeMap must be the same, and, if the matrix is square,
  // getRangeMap and getDomainMap must coincide too.
  //
  const auto nonoverlap_map = this->complex_map();
#ifndef NDEBUG
  TEUCHOS_ASSERT(nonoverlap_map);
#endif
  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  auto graph = Teuchos::rcp(new Tpetra::CrsGraph<int, int>(
      Teuchos::rcp(nonoverlap_map),
      Teuchos::rcp(nonoverlap_map),
      0
      ));

  std::cout << "hasColMap " << graph->hasColMap() << std::endl;
  std::cout << "isLocallyIndexed " << graph->isLocallyIndexed() << std::endl;
  std::cout << "isGloballyIndexed " << graph->isGloballyIndexed() << std::endl;

  const std::vector<edge> edges = this->my_edges();

  // Loop over all edges and put entries wherever two nodes are connected.
  // TODO check if we can use LIDs here
  for (size_t k = 0; k < edges.size(); k++) {
    const Teuchos::Tuple<int,4> & idx = this->edge_lids_complex[k];
    for (int i = 0; i < 4; i++) {
      graph->insertLocalIndices(idx[i], idx);
    }
  }

  graph->fillComplete();

  std::cout << "   build_complex_graph >>" << std::endl;
  return graph;
}
// =============================================================================
Eigen::Vector3d
mesh::
compute_triangle_circumcenter_(
  const std::vector<Eigen::Vector3d> &nodes
  ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(nodes.size(), 3);
#endif
  return this->compute_triangle_circumcenter_(nodes[0], nodes[1], nodes[2]);
}
// =============================================================================
Eigen::Vector3d
mesh::
compute_triangle_circumcenter_(
    const Eigen::Vector3d &node0,
    const Eigen::Vector3d &node1,
    const Eigen::Vector3d &node2
    ) const
{
  Eigen::Vector3d a = node0 - node1;
  Eigen::Vector3d b = node1 - node2;
  Eigen::Vector3d c = node2 - node0;

  const double omega = 2.0 * (a.cross(b)).squaredNorm();

  // don't divide by 0
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      fabs(omega) < 1.0e-10,
      "It seems that the nodes \n"
      << "\n"
      << "   " << node0 << "\n"
      << "   " << node1 << "\n"
      << "   " << node2 << "\n"
      << "\n"
      << "do not form a proper triangle. Abort."
      << std::endl
      );

  const double alpha = - b.dot(b) * a.dot(c) / omega;
  const double beta  = - c.dot(c) * b.dot(a) / omega;
  const double gamma = - a.dot(a) * c.dot(b) / omega;

  return alpha * node0 + beta * node1 + gamma * node2;
}
// =============================================================================
}  // namespace nosh
