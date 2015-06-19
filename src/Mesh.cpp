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
#include "Mesh.hpp"

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
//#include <stk_mesh/base/SkinMesh.hpp>
//
//#include <stk_io/IossBridge.hpp>
//#include <Ioss_SubSystem.h>
////#include <stk_io/MeshReadWriteUtils.hpp>
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
namespace Nosh
{
// =============================================================================
Mesh::
Mesh(
    const std::shared_ptr<const Teuchos::Comm<int>> & _comm,
    const std::shared_ptr<stk::io::StkMeshIoBroker> & broker
    ) :
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  writeTime_(Teuchos::TimeMonitor::getNewTimer("Nosh: Mesh::write")),
  getComplexTime_(Teuchos::TimeMonitor::getNewTimer("Nosh: Mesh::getComplexVector")),
  getMultiTime_(Teuchos::TimeMonitor::getNewTimer("Nosh: Mesh::getMultiVector")),
#endif
  comm(_comm),
  ioBroker_(broker),
  ownedNodes_(this->buildOwnedNodes_(ioBroker_->bulk_data())),
  nodesMap_(this->createEntitiesMap_(ownedNodes_)),
  nodesOverlapMap_(this->createEntitiesMap_(this->getOverlapNodes())),
  complexMap_(this->createComplexMap_(ownedNodes_)),
  complexOverlapMap_(this->createComplexMap_(this->getOverlapNodes())),
  edgeData_(this->createEdgeData_()),
  outputChannel_(0),
  // TODO get index right
  //time_(ioBroker_->get_input_io_region()->get_state_time(index+1)),
  time_(ioBroker_->get_input_io_region()->get_state_time(1)),
  edgeGids(buildEdgeGids_()),
  edgeGidsComplex(buildEdgeGidsComplex_())
{
#ifndef NDEBUG
  // Assert that all processes own nodes
  TEUCHOS_ASSERT_INEQUALITY(ownedNodes_.size(), >, 0);
#endif
}
// =============================================================================
Mesh::
~Mesh()
{
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
Mesh::
complexfield2vector_(
    const ScalarFieldType &realField,
    const ScalarFieldType &imagField
    ) const
{
  // Psi needs to have unique node IDs to be able to compute Norm2().
  // This is required in Belos.
  const auto & ownedNodes = this->getOwnedNodes();

  // Create vector with this respective map.
  auto vector = std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(this->getMapComplex())
      );

  auto vData = vector->getDataNonConst();

#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(vData.size(), 2*ownedNodes.size());
#endif

  // Fill the vector with data from the file.
  for (size_t k = 0; k < ownedNodes.size(); k++) {
    // real part
    double* realVal = stk::mesh::field_data(realField, ownedNodes[k]);
    vData[2*k] = realVal[0];

    // imaginary part
    double* imagVal = stk::mesh::field_data(imagField, ownedNodes[k]);
    vData[2*k+1] = imagVal[0];
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
Mesh::
field2vector_(const ScalarFieldType &field) const
{
  // Get overlap nodes.
  const auto & overlapNodes = this->getOverlapNodes();

  // Create vector with this respective map.
  auto vector = std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(this->getOverlapMap())
      );

  auto vData = vector->getDataNonConst();

#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(vData.size(), 2*overlapNodes.size());
#endif

  // Fill the vector with data from the file.
  for (unsigned int k = 0; k < overlapNodes.size(); k++) {
    double* vals = stk::mesh::field_data(field, overlapNodes[k]);
    // Check if the field is actually there.
#ifndef NDEBUG
    TEUCHOS_ASSERT(vals != NULL);
    //*out << "WARNING: Value for node " << k << " not found.\n" <<
    //  "Probably there is no field given with the state. Using default."
    //  << std::endl;
#endif
    vData[k] = vals[0];
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
Mesh::
field2vector_(
    const VectorFieldType &field,
    const int numComponents
    ) const
{
  // Get overlap nodes.
  const auto & overlapNodes = this->getOverlapNodes();

  // Create vector with this respective map.
  auto vector = std::make_shared<Tpetra::MultiVector<double,int,int>>(
      Teuchos::rcp(this->getOverlapMap()),
      numComponents
      );

  std::vector<Teuchos::ArrayRCP<double>> data(numComponents);
  for (int i = 0; i < numComponents; i++) {
    data[i] = vector->getDataNonConst(i);
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(data[i].size(), overlapNodes.size());
#endif
  }

  // Fill the vector with data from the file.
  for (unsigned int k = 0; k < overlapNodes.size(); k++) {
    const double * const vals =
      stk::mesh::field_data(field, overlapNodes[k]);
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
    for (int i = 0; i < numComponents; i++) {
      data[i][k] = vals[i];
    }
  }

#ifndef NDEBUG
  // Check for NaNs and uninitialized data.
  std::vector<double> r(numComponents);
  // Use NormInf as it's robust against overlapping maps.
  vector->normInf(Teuchos::ArrayView<double>(r));
  for (int i = 0; i < numComponents; i++) {
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
Mesh::
openFile(
    const std::string &outputFile
    )
{
  outputChannel_ = ioBroker_->create_output_mesh(
      outputFile,
      stk::io::WRITE_RESULTS
      );
  const stk::mesh::FieldVector &fields = ioBroker_->meta_data().get_fields();
  for (size_t i=0; i < fields.size(); i++) {
    if (*stk::io::get_field_role(*fields[i]) == Ioss::Field::TRANSIENT) {
      ioBroker_->add_field(outputChannel_, *fields[i]);
    }
  }

  return;
}
// =============================================================================
void
Mesh::
write(const double time) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*writeTime_);
#endif

  // Write it out to the file that's been specified in mesh_.
  // The methods returns the output step (but we ignore it).
  //(void) ioBroker_->process_output_request(
  //    outputChannel_,
  //    time
  //    );
  static int step = 0;
  (void) ioBroker_->process_output_request(outputChannel_, step++);

  return;
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
Mesh::
getVector(const std::string & fieldName) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(ioBroker_);
#endif
  const ScalarFieldType * const field =
    ioBroker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        fieldName
        );

  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      field == NULL,
      "Scalar field \"" << fieldName << "\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  return this->field2vector_(*field);
}
// =============================================================================
std::shared_ptr<Tpetra::MultiVector<double,int,int>>
Mesh::
getMultiVector(const std::string & fieldName) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*getMultiTime_);
#endif
#ifndef NDEBUG
  TEUCHOS_ASSERT(ioBroker_);
#endif

  const VectorFieldType * const field =
    ioBroker_->bulk_data().mesh_meta_data().get_field<VectorFieldType>(
        stk::topology::NODE_RANK,
        fieldName
        );

  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      field == NULL,
      "Vector field \"" << fieldName << "\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  // TODO remove the hardcoded "3"
  return this->field2vector_(*field, 3);
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
Mesh::
getComplexVector(const std::string & fieldName) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*getComplexTime_);
#endif
#ifndef NDEBUG
  TEUCHOS_ASSERT(ioBroker_);
#endif
  const ScalarFieldType * const r_field =
    ioBroker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        fieldName + "_R"
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      r_field == NULL,
      "Scalar field \"" << fieldName << "_R\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  const ScalarFieldType * const i_field =
    ioBroker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        fieldName + "_Z"
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      i_field == NULL,
      "Scalar field \"" << fieldName << "_Z\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  return this->complexfield2vector_(*r_field, *i_field);
}
// =============================================================================
void
Mesh::
insertVector(
    const Tpetra::Vector<double,int,int> & x,
    const std::string & fieldName
    ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*writeTime_);
#endif

#ifndef NDEBUG
  TEUCHOS_ASSERT(ioBroker_);
#endif

  ScalarFieldType * x_field =
    ioBroker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        fieldName
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      x_field == NULL,
      "Scalar field \"" << fieldName << "\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  // Zero out all nodal values, including the overlaps.
  const auto & overlapNodes = this->getOverlapNodes();
  for (unsigned int k = 0; k < overlapNodes.size(); k++) {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data(*x_field, overlapNodes[k]);
    *localPsiR = 0.0;
  }

  auto xData = x.getData();

  // Set owned nodes.
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(xData.size(), ownedNodes_.size());
#endif
  for (unsigned int k = 0; k < ownedNodes_.size(); k++) {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data(*x_field, ownedNodes_[k]);
    *localPsiR = xData[k];
  }

  // This communication updates the field values on un-owned nodes it is
  // correct because the zeroSolutionField above zeros them all and the
  // getSolutionField only sets the owned nodes.
  // TODO combine these fields into a vector of fields
  std::vector<stk::mesh::FieldBase*> tmp(1, x_field);
  stk::mesh::parallel_sum(ioBroker_->bulk_data(), tmp);

  return;
}
// =============================================================================
void
Mesh::
insertComplexVector(
    const Tpetra::Vector<double,int,int> & psi,
    const std::string & fieldName
    ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*writeTime_);
#endif

#ifndef NDEBUG
  TEUCHOS_ASSERT(ioBroker_);
#endif
  ScalarFieldType * psir_field =
    ioBroker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        fieldName + "_R"
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      psir_field == NULL,
      "Scalar field \"" << fieldName << "_R\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  ScalarFieldType * psii_field =
    ioBroker_->bulk_data().mesh_meta_data().get_field<ScalarFieldType>(
        stk::topology::NODE_RANK,
        fieldName + "_Z"
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      psii_field == NULL,
      "Scalar field \"" << fieldName << "_Z\" not found in database. "
      << "Is it present in the input file at all? Check with io_info."
      );

  // Zero out all nodal values, including the overlaps.
  const auto & overlapNodes = this->getOverlapNodes();
  for (unsigned int k = 0; k < overlapNodes.size(); k++) {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data(*psir_field, overlapNodes[k]);
    *localPsiR = 0.0;
    double* localPsiI = stk::mesh::field_data(*psii_field, overlapNodes[k]);
    *localPsiI = 0.0;
  }

  auto psiData = psi.getData();

  // Set owned nodes.
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(psiData.size(), 2*ownedNodes_.size());
#endif
  for (unsigned int k = 0; k < ownedNodes_.size(); k++) {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data(*psir_field, ownedNodes_[k]);
    *localPsiR = psiData[2*k];
    double* localPsiI = stk::mesh::field_data(*psii_field, ownedNodes_[k]);
    *localPsiI = psiData[2*k+1];
  }

  // This communication updates the field values on un-owned nodes
  // it is correct because the zeroSolutionField above zeros them all
  // and the getSolutionField only sets the owned nodes.
  // TODO combine these fields into a vector of fields
  std::vector<stk::mesh::FieldBase*> tmp(1, psir_field);
  stk::mesh::parallel_sum(ioBroker_->bulk_data(), tmp);
  std::vector<stk::mesh::FieldBase*> tmp2(1, psii_field);
  stk::mesh::parallel_sum(ioBroker_->bulk_data(), tmp2);

  return;
}

// =============================================================================
std::vector<stk::mesh::Entity>
Mesh::
getOwnedCells() const
{
  // get owned elements
  stk::mesh::Selector select_owned_in_part =
    stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().universal_part())
    & stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().locally_owned_part());
  std::vector<stk::mesh::Entity> cells;
  stk::mesh::get_selected_entities(
      select_owned_in_part,
      ioBroker_->bulk_data().buckets(stk::topology::ELEMENT_RANK),
      cells
      );
  return cells;
}
// =============================================================================
std::vector<stk::mesh::Entity>
Mesh::
getOverlapEdges() const
{
  // get overlap edges
  stk::mesh::Selector select_overlap_in_part =
    stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().universal_part())
    & (stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().locally_owned_part())
       |stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().globally_shared_part()));

  std::vector<stk::mesh::Entity> edges;
  stk::mesh::get_selected_entities(
      select_overlap_in_part,
      ioBroker_->bulk_data().buckets(stk::topology::EDGE_RANK),
      edges
      );
  return edges;
}
// =============================================================================
const VectorFieldType &
Mesh::
getNodeField(const std::string & fieldName) const {
  const VectorFieldType * const field =
    ioBroker_->bulk_data().mesh_meta_data().get_field<VectorFieldType>(
        stk::topology::NODE_RANK,
        fieldName
        );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      field == NULL,
      "Invalid field name \"" << fieldName << "\"."
      );
  return *field;
}
// =============================================================================
std::vector<stk::mesh::Entity>
Mesh::
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
Mesh::
getOverlapNodes() const
{
  //  overlapnodes used for overlap map -- stored for changing coords
  stk::mesh::Selector select_overlap_in_part =
    stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().universal_part())
    & (stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().locally_owned_part())
       |stk::mesh::Selector(ioBroker_->bulk_data().mesh_meta_data().globally_shared_part()));

  std::vector<stk::mesh::Entity> overlapNodes;
  stk::mesh::get_selected_entities(
      select_overlap_in_part,
      ioBroker_->bulk_data().buckets(stk::topology::NODE_RANK),
      overlapNodes
      );

  return overlapNodes;
}
// =============================================================================
const std::vector<Teuchos::Tuple<int,2>>
Mesh::
buildEdgeGids_() const
{
  const std::vector<edge> edges = this->getMyEdges();

  std::vector<Teuchos::Tuple<int,2>> gic(edges.size());

  int gidT0, gidT1;
  for (std::size_t k = 0; k < edges.size(); k++) {
    gidT0 = this->gid(std::get<0>(edges[k]));
    gidT1 = this->gid(std::get<1>(edges[k]));
    gic[k] = Teuchos::tuple(gidT0, gidT1);
  }

  return gic;
}
// =============================================================================
const std::vector<Teuchos::Tuple<int,4>>
Mesh::
buildEdgeGidsComplex_() const
{
  const std::vector<edge> edges = this->getMyEdges();

  std::vector<Teuchos::Tuple<int,4>> gic(edges.size());

  int gidT0, gidT1;
  for (std::size_t k = 0; k < edges.size(); k++) {
    gidT0 = this->gid(std::get<0>(edges[k]));
    gidT1 = this->gid(std::get<1>(edges[k]));
    gic[k] = Teuchos::tuple(2*gidT0, 2*gidT0 + 1, 2*gidT1, 2*gidT1 + 1);
  }

  return gic;
}
// =============================================================================
std::shared_ptr<const Tpetra::Map<int,int>>
Mesh::
createEntitiesMap_(const std::vector<stk::mesh::Entity> &entityList) const
{
  const size_t numEntities = entityList.size();
  std::vector<int> gids(numEntities);
  for (size_t i = 0; i < numEntities; i++) {
    gids[i] = ioBroker_->bulk_data().identifier(entityList[i]) - 1;
  }

  return std::make_shared<Tpetra::Map<int,int>>(
      -1, Teuchos::ArrayView<int>(gids), 0, Teuchos::rcp(this->comm)
      );
}
// =============================================================================
std::shared_ptr<const Tpetra::Map<int,int>>
Mesh::
createComplexMap_(const std::vector<stk::mesh::Entity> &nodeList) const
{
  // Create a map for real/imaginary out of this.
  const size_t numDof = 2 * nodeList.size();
  std::vector<int> gids(numDof);
  for (size_t k = 0; k < nodeList.size(); k++) {
    int globalNodeId = ioBroker_->bulk_data().identifier(nodeList[k]) - 1;
    gids[2*k]   = 2*globalNodeId;
    gids[2*k+1] = 2*globalNodeId + 1;
  }

  return std::make_shared<Tpetra::Map<int,int>>(
      -1, Teuchos::ArrayView<int>(gids), 0, Teuchos::rcp(this->comm)
      );
}
// =============================================================================
Mesh::EdgesContainer
Mesh::
createEdgeData_()
{
  std::vector<stk::mesh::Entity> cells = this->getOwnedCells();
  size_t numLocalCells = cells.size();

  Mesh::EdgesContainer edgeData = {
    // Local edge ID -> Local node IDs.
    std::vector<std::tuple<stk::mesh::Entity, stk::mesh::Entity> >(),
    // Local cell ID -> Local edge IDs.
    std::vector<std::vector<int>>(numLocalCells)
    };

  // Used to determine if an edge sits on the domain boundary
  std::vector<size_t> numCellsPerEdge();

  // This std::map keeps track of how nodes and edges are connected. If
  // `nodeEdges((3,4)) == 17`, then the nodes (3,4) are connected by edge 17.
  // Unfortunately, std::tuples can't be compared with '<'. Provide a function
  // pointer that implements lexicographic comparison.
  // See http://www.cplusplus.com/reference/stl/map/map/.
  std::map<std::tuple<stk::mesh::Entity, stk::mesh::Entity>, int> nodesEdge;

  //const EntityComp ec(ioBroker_->bulk_data());

  // Loop over all owned cells.
  unsigned int edgeLID = 0;
  for (size_t cellLID = 0; cellLID < numLocalCells; cellLID++) {
    // Loop over all pairs of local nodes.
    stk::mesh::Entity const * localNodes
      = ioBroker_->bulk_data().begin_nodes(cells[cellLID]);
    const size_t numLocalNodes =
      ioBroker_->bulk_data().num_nodes(cells[cellLID]);

    //stk::mesh::PairIterRelation nodesIterator =
    //  cells[cellLID]->relations(metaData.node_rank());
    //unsigned int numLocalNodes = nodesIterator.size();
    size_t numLocalEdges = numLocalNodes*(numLocalNodes-1) / 2;

    edgeData.cellEdges[cellLID] = std::vector<int>(numLocalEdges);

    // Gather the node entities.
    std::vector<stk::mesh::Entity> nodes(numLocalNodes);
    for (size_t k = 0; k < numLocalNodes; k++) {
      nodes[k] = localNodes[k];
    }

    // Sort nodes. This is necessary to make sure that the tuples formed below
    // are always sorted such they are unique keys (and {3,7}, {7,3} are
    // recognized as the same edge).
    std::sort(nodes.begin(), nodes.end());

    // In a simplex, the edges are exactly the connection between each pair of
    // nodes. Hence, loop over pairs of nodes.
    unsigned int edgeIndex = 0;
    edge edgeNodes;
    for (size_t e0 = 0; e0 < numLocalNodes; e0++) {
      std::get<0>(edgeNodes) = nodes[e0];
      for (size_t e1 = e0+1; e1 < numLocalNodes; e1++) {
        std::get<1>(edgeNodes) = nodes[e1];
        // As nodes are sorted and by their identifiers, edgeNodes are sorted
        // too. This is necessary as otherwise the edge {3,7} could not be
        // identified as {7,3}.
        // Check if edgeNodes is in the map.
        auto it = nodesEdge.find(edgeNodes);
        if (it != nodesEdge.end()) {
          // Edge is already accounted for.
          edgeData.cellEdges[cellLID][edgeIndex] = it->second;
        } else {  // Edge not found -- insert it.
          nodesEdge[edgeNodes] = edgeLID; // for householding in this method
          edgeData.edgeNodes.push_back(edgeNodes); // for looping over edges
          edgeData.cellEdges[cellLID][edgeIndex] = edgeLID; // for this->computeEdgeCoefficients_
          edgeLID++;
        }
        edgeIndex++;
      }
    }
  }

  return edgeData;
}
// =============================================================================
Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
Mesh::
buildGraph() const
{
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from
  // the ownedNodes map, and the overlapping one from the
  // overlapNodes.
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
  const auto noMap = this->getMap();
#ifndef NDEBUG
  TEUCHOS_ASSERT(noMap);
#endif
  auto graph = Tpetra::createCrsGraph(Teuchos::rcp(noMap));

  const std::vector<edge> edges = this->getMyEdges();

  // Loop over all edges and put entries wherever two nodes are connected.
  // TODO check if we can use LIDs here
  for (size_t k = 0; k < edges.size(); k++) {
    const Teuchos::Tuple<int,2> & idx = this->edgeGids[k];
    for (int i = 0; i < 2; i++) {
      graph->insertGlobalIndices(idx[i], idx);
    }
  }

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  // TODO specify noMap?
  graph->fillComplete();

  return graph;
}
// =============================================================================
Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
Mesh::
buildComplexGraph() const
{
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from
  // the ownedNodes map, and the overlapping one from the
  // overlapNodes.
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
  const auto noMap = this->getMapComplex();
#ifndef NDEBUG
  TEUCHOS_ASSERT(noMap);
#endif
  auto graph = Tpetra::createCrsGraph(Teuchos::rcp(noMap));

  const std::vector<edge> edges = this->getMyEdges();

  // Loop over all edges and put entries wherever two nodes are connected.
  // TODO check if we can use LIDs here
  for (size_t k = 0; k < edges.size(); k++) {
    const Teuchos::Tuple<int,4> & idx = this->edgeGidsComplex[k];
    for (int i = 0; i < 4; i++) {
      graph->insertGlobalIndices(idx[i], idx);
    }
  }

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  // TODO specify noMap?
  graph->fillComplete();

  return graph;
}
// =============================================================================
Eigen::Vector3d
Mesh::
computeTriangleCircumcenter_(
  const std::vector<Eigen::Vector3d> &nodes
  ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(nodes.size(), 3);
#endif
  return this->computeTriangleCircumcenter_(nodes[0], nodes[1], nodes[2]);
}
// =============================================================================
Eigen::Vector3d
Mesh::
computeTriangleCircumcenter_(
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
}  // namespace Nosh
