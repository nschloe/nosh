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

//#include <map>
//#include <string>
//#include <algorithm>
//#include <vector>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <stk_mesh/base/MetaData.hpp>
//#include <stk_mesh/base/BulkData.hpp>
//#include <stk_mesh/base/Entity.hpp>
//// #include <stk_mesh/base/Field.hpp>
//#include <stk_mesh/base/FieldBase.hpp>
//// #include <stk_mesh/base/Comm.hpp> // for comm_mesh_counts
//#include <stk_mesh/base/CreateAdjacentEntities.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldParallel.hpp> // for parallel_sum
//#include <stk_mesh/base/MetaData.hpp>
//#include <stk_mesh/base/Comm.hpp>

//#include <stk_mesh/base/CreateEdges.hpp>
//#include <stk_mesh/base/CreateFaces.hpp>
//#include <stk_mesh/base/Skinmesh.hpp>
//
//#include <stk_io/IossBridge.hpp>
//#include <Ioss_SubSystem.h>
////#include <stk_io/mesh_readWriteUtils.hpp>
//#include <Ionit_Initializer.h>
//#include <Ioss_IOFactory.h>
#include <Ioss_Region.h>

//#ifdef HAVE_MPI
//// Rebalance
//#include <stk_rebalance/Rebalance.hpp>
//#include <stk_rebalance_utils/RebalanceUtils.hpp>
//#include <stk_rebalance/Partition.hpp>
//#include <stk_rebalance/ZoltanPartition.hpp>
//#endif

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif
namespace nosh
{
// =============================================================================
mesh::
mesh(
    const std::shared_ptr<const Teuchos::Comm<int>> & _comm,
    const std::shared_ptr<stk::io::StkMeshIoBroker> & broker
    ) :
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  write_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: mesh::write")),
  getComplex_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: mesh::get_complex_vector")),
  multi_time_(Teuchos::TimeMonitor::getNewTimer("Nosh: mesh::get_multi_vector")),
#endif
  comm(_comm),
  io_broker_(broker),
  owned_nodes_(this->buildOwnedNodes_(io_broker_->bulk_data())),
  nodes_map_(this->buildEntitiesMap_(owned_nodes_)),
  nodes_overlap_map_(this->buildEntitiesMap_(this->get_overlap_nodes())),
  complex_map_(this->buildComplexMap_(owned_nodes_)),
  complex_overlap_map_(this->buildComplexMap_(this->get_overlap_nodes())),
  edge_data_(this->buildEdge_data_()),
  output_channel_(0),
  // TODO get index right
  //time_(io_broker_->get_input_io_region()->get_state_time(index+1)),
  time_(io_broker_->get_input_io_region()->get_state_time(1)),
  edge_gids(buildEdgeGids_()),
  edge_gids_complex(buildEdgeGidsComplex_())
{
#ifndef NDEBUG
  // Assert that all processes own nodes
  TEUCHOS_ASSERT_INEQUALITY(owned_nodes_.size(), >, 0);
#endif
}
// =============================================================================
mesh::
~mesh()
{
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh::
complexfield_to_vector_(
    const ScalarFieldType &real_field,
    const ScalarFieldType &imag_field
    ) const
{
  // Psi needs to have unique node IDs to be able to compute Norm2().
  // This is required in Belos.
  const auto & _owned_nodes = this->owned_nodes();

  // Create vector with this respective map.
  auto vector = std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(this->complex_map())
      );

  auto v_data = vector->getDataNonConst();

#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(v_data.size(), 2 * _owned_nodes.size());
#endif

  // Fill the vector with data from the file.
  for (size_t k = 0; k < _owned_nodes.size(); k++) {
    // real part
    double* real_val = stk::mesh::field_data(real_field, _owned_nodes[k]);
    v_data[2*k] = real_val[0];

    // imaginary part
    double* imag_val = stk::mesh::field_data(imag_field, _owned_nodes[k]);
    v_data[2*k+1] = imag_val[0];
  }

#ifndef NDEBUG
  // Use NormInf as it's robust against overlapping maps.
  const double r = vector->normInf();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      r != r || r > 1.0e100,
      "The input data seems flawed. Abort."
      );
#endif

  return vector;
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh::
field_to_vector_(const ScalarFieldType &field) const
{
  // Get overlap nodes.
  const auto & overlap_nodes = this->get_overlap_nodes();

  // Create vector with this respective map.
  auto vector = std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(this->overlap_map())
      );

  auto v_data = vector->getDataNonConst();

#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(v_data.size(), 2*overlap_nodes.size());
#endif

  // Fill the vector with data from the file.
  for (unsigned int k = 0; k < overlap_nodes.size(); k++) {
    double* vals = stk::mesh::field_data(field, overlap_nodes[k]);
    // Check if the field is actually there.
#ifndef NDEBUG
    TEUCHOS_ASSERT(vals != NULL);
    //*out << "WARNING: _value for node " << k << " not found.\n" <<
    //  "Probably there is no field given with the state. Using default."
    //  << std::endl;
#endif
    v_data[k] = vals[0];
  }

#ifndef NDEBUG
  // Use NormInf as it's robust against overlapping maps.
  const double r = vector->normInf();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      r != r || r > 1.0e100,
      "The input data seems flawed. Abort."
      );
#endif

  return vector;
}
// =============================================================================
std::shared_ptr<Tpetra::MultiVector<double,int,int>>
mesh::
field_to_vector_(
    const vector_fieldType &field,
    const int num_components
    ) const
{
  // Get overlap nodes.
  const auto & overlap_nodes = this->get_overlap_nodes();

  // Create vector with this respective map.
  auto vector = std::make_shared<Tpetra::MultiVector<double,int,int>>(
      Teuchos::rcp(this->overlap_map()),
      num_components
      );

  std::vector<Teuchos::ArrayRCP<double>> data(num_components);
  for (int i = 0; i < num_components; i++) {
    data[i] = vector->getDataNonConst(i);
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(data[i].size(), overlap_nodes.size());
#endif
  }

  // Fill the vector with data from the file.
  for (unsigned int k = 0; k < overlap_nodes.size(); k++) {
    const double * const vals =
      stk::mesh::field_data(field, overlap_nodes[k]);
#ifndef NDEBUG
    // Check if the field is actually there.
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        vals == NULL,
        "Field value for node " << k << " not found.\n" <<
        "Probably there is no field given with the state."
        );
#endif
    // Copy over.
    // A multivector isn't actually a good data structure for this.  What would
    // be needed is a vector where each entry has k components. This way, the
    // data could stick together.
    for (int i = 0; i < num_components; i++) {
      data[i][k] = vals[i];
    }
  }

#ifndef NDEBUG
  // Check for NaNs and uninitialized data.
  std::vector<double> r(num_components);
  // Use NormInf as it's robust against overlapping maps.
  vector->normInf(Teuchos::ArrayView<double>(r));
  for (int i = 0; i < num_components; i++) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        r[i] != r[i] || r[i] > 1.0e100,
        "The input data seems flawed. Abort."
        );
  }
#endif

  return vector;
}
// =============================================================================
void
mesh::
open_file(
    const std::string &output_file
    )
{
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

  return;
}
// =============================================================================
void
mesh::
write(const double _time) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*write_time_);
#endif

  (void) _time;
  // Write it out to the file that's been specified in mesh_.
  // The methods returns the output step (but we ignore it).
  //(void) io_broker_->process_output_request(
  //    output_channel_,
  //    _time
  //    );
  static int step = 0;
  (void) io_broker_->process_output_request(output_channel_, step++);

  return;
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh::
get_vector(const std::string & field_name) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(io_broker_);
#endif
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
}
// =============================================================================
std::shared_ptr<Tpetra::MultiVector<double,int,int>>
mesh::
get_multi_vector(const std::string & field_name) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*multi_time_);
#endif
#ifndef NDEBUG
  TEUCHOS_ASSERT(io_broker_);
#endif

  const vector_fieldType * const field =
    io_broker_->bulk_data().mesh_meta_data().get_field<vector_fieldType>(
        stk::topology::NODE_RANK,
        field_name
        );

  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      field == NULL,
      "Vector field \"" << field_name << "\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  // TODO remove the hardcoded "3"
  return this->field_to_vector_(*field, 3);
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh::
get_complex_vector(const std::string & field_name) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*getComplex_time_);
#endif
#ifndef NDEBUG
  TEUCHOS_ASSERT(io_broker_);
#endif
  const ScalarFieldType * const r_field =
    io_broker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        field_name + "_R"
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      r_field == NULL,
      "Scalar field \"" << field_name << "_R\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  const ScalarFieldType * const i_field =
    io_broker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        field_name + "_Z"
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      i_field == NULL,
      "Scalar field \"" << field_name << "_Z\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  return this->complexfield_to_vector_(*r_field, *i_field);
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

#ifndef NDEBUG
  TEUCHOS_ASSERT(io_broker_);
#endif

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

  return;
}
// =============================================================================
void
mesh::
insert_complex_vector(
    const Tpetra::Vector<double,int,int> & psi,
    const std::string & field_name
    ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*write_time_);
#endif

#ifndef NDEBUG
  TEUCHOS_ASSERT(io_broker_);
#endif
  ScalarFieldType * psi_r_field =
    io_broker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        field_name + "_R"
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      psi_r_field == NULL,
      "Scalar field \"" << field_name << "_R\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  ScalarFieldType * psi_i_field =
    io_broker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        field_name + "_Z"
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      psi_i_field == NULL,
      "Scalar field \"" << field_name << "_Z\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  // Zero out all nodal values, including the overlaps.
  const auto & overlap_nodes = this->get_overlap_nodes();
  for (unsigned int k = 0; k < overlap_nodes.size(); k++) {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data(*psi_r_field, overlap_nodes[k]);
    *localPsiR = 0.0;
    double* localPsiI = stk::mesh::field_data(*psi_i_field, overlap_nodes[k]);
    *localPsiI = 0.0;
  }

  auto psi_data = psi.getData();

  // Set owned nodes.
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(psi_data.size(), 2*owned_nodes_.size());
#endif
  for (unsigned int k = 0; k < owned_nodes_.size(); k++) {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data(*psi_r_field, owned_nodes_[k]);
    *localPsiR = psi_data[2*k];
    double* localPsiI = stk::mesh::field_data(*psi_i_field, owned_nodes_[k]);
    *localPsiI = psi_data[2*k+1];
  }

  // This communication updates the field values on un-owned nodes
  // it is correct because the zeroSolutionField above zeros them all
  // and the getSolutionField only sets the owned nodes.
  // TODO combine these fields into a vector of fields
  std::vector<stk::mesh::FieldBase*> tmp(1, psi_r_field);
  stk::mesh::parallel_sum(io_broker_->bulk_data(), tmp);
  std::vector<stk::mesh::FieldBase*> tmp2(1, psi_i_field);
  stk::mesh::parallel_sum(io_broker_->bulk_data(), tmp2);

  return;
}

// =============================================================================
std::vector<stk::mesh::Entity>
mesh::
get_owned_cells() const
{
  // get owned elements
  stk::mesh::Selector select_owned_in_part =
    stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().universal_part())
    & stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().locally_owned_part());
  std::vector<stk::mesh::Entity> cells;
  stk::mesh::get_selected_entities(
      select_owned_in_part,
      io_broker_->bulk_data().buckets(stk::topology::ELEMENT_RANK),
      cells
      );
  return cells;
}
// =============================================================================
std::vector<stk::mesh::Entity>
mesh::
get_overlap_edges() const
{
  // get overlap edges
  stk::mesh::Selector select_overlap_in_part =
    stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().universal_part())
    & (stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().locally_owned_part())
       |stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().globally_shared_part()));

  std::vector<stk::mesh::Entity> edges;
  stk::mesh::get_selected_entities(
      select_overlap_in_part,
      io_broker_->bulk_data().buckets(stk::topology::EDGE_RANK),
      edges
      );
  return edges;
}
// =============================================================================
const vector_fieldType &
mesh::
get_node_field(const std::string & field_name) const {
  const vector_fieldType * const field =
    io_broker_->bulk_data().mesh_meta_data().get_field<vector_fieldType>(
        stk::topology::NODE_RANK,
        field_name
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      field == NULL,
      "Invalid field name \"" << field_name << "\"."
      );
  return *field;
}
// =============================================================================
std::vector<stk::mesh::Entity>
mesh::
buildOwnedNodes_(const stk::mesh::BulkData & myBulkData) const
{
  stk::mesh::Selector select_owned_in_part =
    stk::mesh::Selector(myBulkData.mesh_meta_data().universal_part())
    & stk::mesh::Selector(myBulkData.mesh_meta_data().locally_owned_part());

  std::vector<stk::mesh::Entity> on;
  stk::mesh::get_selected_entities(
      select_owned_in_part,
      myBulkData.buckets(stk::topology::NODE_RANK),
      on
      );
  return on;
}
// =============================================================================
std::vector<stk::mesh::Entity>
mesh::
get_overlap_nodes() const
{
  //  overlapnodes used for overlap map -- stored for changing coords
  stk::mesh::Selector select_overlap_in_part =
    stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().universal_part())
    & (stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().locally_owned_part())
       |stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().globally_shared_part()));

  std::vector<stk::mesh::Entity> overlap_nodes;
  stk::mesh::get_selected_entities(
      select_overlap_in_part,
      io_broker_->bulk_data().buckets(stk::topology::NODE_RANK),
      overlap_nodes
      );

  return overlap_nodes;
}
// =============================================================================
const std::vector<Teuchos::Tuple<int,2>>
mesh::
buildEdgeGids_() const
{
  const std::vector<edge> edges = this->my_edges();

  std::vector<Teuchos::Tuple<int,2>> _edge_gids(edges.size());

  int gidT0, gidT1;
  for (std::size_t k = 0; k < edges.size(); k++) {
    gidT0 = this->gid(std::get<0>(edges[k]));
    gidT1 = this->gid(std::get<1>(edges[k]));
    _edge_gids[k] = Teuchos::tuple(gidT0, gidT1);
  }

  return _edge_gids;
}
// =============================================================================
const std::vector<Teuchos::Tuple<int,4>>
mesh::
buildEdgeGidsComplex_() const
{
  const std::vector<edge> edges = this->my_edges();

  std::vector<Teuchos::Tuple<int,4>> _edge_gids_complex(edges.size());

  int gidT0, gidT1;
  for (std::size_t k = 0; k < edges.size(); k++) {
    gidT0 = this->gid(std::get<0>(edges[k]));
    gidT1 = this->gid(std::get<1>(edges[k]));
    _edge_gids_complex[k] =
      Teuchos::tuple(2*gidT0, 2*gidT0 + 1, 2*gidT1, 2*gidT1 + 1);
  }

  return _edge_gids_complex;
}
// =============================================================================
std::shared_ptr<const Tpetra::Map<int,int>>
mesh::
buildEntitiesMap_(const std::vector<stk::mesh::Entity> &entityList) const
{
  const size_t numEntities = entityList.size();
#ifndef NDEBUG
  TEUCHOS_ASSERT_INEQUALITY(numEntities, >, 0);
#endif
  std::vector<int> gids(numEntities);
  for (size_t i = 0; i < numEntities; i++) {
    gids[i] = io_broker_->bulk_data().identifier(entityList[i]) - 1;
  }

  return std::make_shared<Tpetra::Map<int,int>>(
      -1, Teuchos::ArrayView<int>(gids), 0, Teuchos::rcp(this->comm)
      );
}
// =============================================================================
std::shared_ptr<const Tpetra::Map<int,int>>
mesh::
buildComplexMap_(const std::vector<stk::mesh::Entity> &node_list) const
{
  // Create a map for real/imaginary out of this.
  const size_t num_dof = 2 * node_list.size();
  std::vector<int> gids(num_dof);
  for (size_t k = 0; k < node_list.size(); k++) {
    int global_node_id = io_broker_->bulk_data().identifier(node_list[k]) - 1;
    gids[2*k]   = 2 * global_node_id;
    gids[2*k+1] = 2 * global_node_id + 1;
  }

  return std::make_shared<Tpetra::Map<int,int>>(
      -1, Teuchos::ArrayView<int>(gids), 0, Teuchos::rcp(this->comm)
      );
}
// =============================================================================
mesh::edges_container
mesh::
buildEdge_data_()
{
  std::vector<stk::mesh::Entity> cells = this->get_owned_cells();
  size_t num_local_cells = cells.size();

  mesh::edges_container edge_data = {
    // Local edge ID -> Local node IDs.
    std::vector<std::tuple<stk::mesh::Entity, stk::mesh::Entity> >(),
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
  std::map<std::tuple<stk::mesh::Entity, stk::mesh::Entity>, int> nodesEdge;

  //const EntityComp ec(io_broker_->bulk_data());

  // Loop over all owned cells.
  unsigned int edge_local_id = 0;
  for (size_t cellLID = 0; cellLID < num_local_cells; cellLID++) {
    // Loop over all pairs of local nodes.
    stk::mesh::Entity const * local_nodes
      = io_broker_->bulk_data().begin_nodes(cells[cellLID]);
    const size_t num_local_nodes =
      io_broker_->bulk_data().num_nodes(cells[cellLID]);

    //stk::mesh::PairIterRelation nodesIterator =
    //  cells[cellLID]->relations(meta_data.node_rank());
    //unsigned int num_local_nodes = nodesIterator.size();
    size_t num_local_edges = num_local_nodes * (num_local_nodes - 1) / 2;

    edge_data.cell_edges[cellLID] = std::vector<int>(num_local_edges);

    // Gather the node entities.
    std::vector<stk::mesh::Entity> nodes(num_local_nodes);
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

  return edge_data;
}
// =============================================================================
Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
mesh::
build_graph() const
{
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
  auto graph = Tpetra::createCrsGraph(Teuchos::rcp(nonoverlap_map));

  const std::vector<edge> edges = this->my_edges();

  // Loop over all edges and put entries wherever two nodes are connected.
  // TODO check if we can use LIDs here
  for (size_t k = 0; k < edges.size(); k++) {
    const Teuchos::Tuple<int,2> & idx = this->edge_gids[k];
    for (int i = 0; i < 2; i++) {
      graph->insertGlobalIndices(idx[i], idx);
    }
  }

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  // TODO specify nonoverlap_map?
  graph->fillComplete();

  return graph;
}
// =============================================================================
Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
mesh::
build_complex_graph() const
{
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from the
  // owned_nodes map, and the overlapping one from the overlap_nodes.
  // Let's illustrate the implications with the example of the matrix
  //   [ 2 1   ]
  //   [ 1 2 1 ]
  //   [   1 2 ].
  // Suppose subdomain 1 consists of node 1, subdomain 2 of node 3, and node 2
  // forms the boundary between them.
  // For two processes, if process 1 owns nodes 1 and 2, the matrix will be
  // split as
  //   [ 2 1   ]   [       ]
  //   [ 1 2 1 ] + [       ]
  //   [       ]   [   1 2 ].
  // The vectors always need to have a unique map (otherwise, norms cannot be
  // computed by Epetra), so let's assume they have the map ([1, 2], [3]).
  // The communucation for a matrix-vector multiplication Ax=y needs to be:
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
  // In the general case, assuming that the number of nodes adjacent to a
  // boundary (on one side) are approximately the number of nodes on that
  // boundary, there is not much difference in communication between the
  // patterns.
  // What does differ, though, is the workload on the processes during the
  // computation phase: Process 1 that owns the whole boundary, has to compute
  // more than process 2.
  // Notice, however, that the total number of computations is lower in
  // scenario 1 (7 vs. 8 FLOPs); the same is true for storage.
  // Hence, it comes down to the question whether or not the mesh generator
  // provided a fair share of the boundary nodes.  If yes, then scenario 1 will
  // yield approximately even computation times; if not, then scenario 2 will
  // guarantee equal computation times at the cost of higher total storage and
  // computation needs.
  //
  // Remark:
  // This matrix will later be fed into ML. ML has certain restrictions as to
  // what maps can be used. One of those is that RowMatrixRowMap() and
  // getRangeMap must be the same, and, if the matrix is square, getRangeMap
  // and getDomainMap must coincide too.
  //
  const auto nonoverlap_map = this->complex_map();
#ifndef NDEBUG
  TEUCHOS_ASSERT(nonoverlap_map);
#endif
  auto graph = Tpetra::createCrsGraph(Teuchos::rcp(nonoverlap_map));

  const std::vector<edge> edges = this->my_edges();

  // Loop over all edges and put entries wherever two nodes are connected.
  // TODO check if we can use LIDs here
  for (size_t k = 0; k < edges.size(); k++) {
    const Teuchos::Tuple<int,4> & idx = this->edge_gids_complex[k];
    for (int i = 0; i < 4; i++) {
      graph->insertGlobalIndices(idx[i], idx);
    }
  }

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  // TODO specify nonoverlap_map?
  graph->fillComplete();

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
