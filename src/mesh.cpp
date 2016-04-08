#include "mesh.hpp"

#include <moab/Core.hpp>
#include <moab/ParallelComm.hpp>
#include <MBParallelConventions.h>
#include <moab/Skinner.hpp>

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
  mbw_(std::make_shared<moab_wrap>(mb)),
  mcomm_(mcomm),
  vertices_map_(this->get_map_(this->get_owned_gids_())),
  vertices_overlap_map_(this->get_map_(this->get_overlap_gids_())),
  complex_map_(this->get_map_(this->complexify_(this->get_owned_gids_()))),
  complex_overlap_map_(
    this->get_map_(this->complexify_(this->get_overlap_gids_()))
    )
  ,relations_(this->build_entity_relations_())
  ,edge_lids(build_edge_lids_())
  ,edge_lids_complex(build_edge_lids_complex_())
  ,edge_gids(build_edge_gids_())
  ,edge_gids_complex(build_edge_gids_complex_())
  ,boundary_skin_(compute_boundary_skin_())
  ,boundary_vertices(compute_boundary_vertices_(boundary_skin_))
  ,meshsets_(create_default_meshsets_())
{
// TODO resurrect
//#ifndef NDEBUG
//  // Assert that all processes own vertices
//  TEUCHOS_ASSERT_INEQUALITY(owned_vertices_.size(), >, 0);
//#endif
}
// =============================================================================
mesh::
~mesh()
{
}
// =============================================================================
std::map<std::string, moab::EntityHandle>
mesh::
create_default_meshsets_()
{
  // create a meshset for all boundary vertices
  const auto boundary = this->mbw_->create_meshset(moab::MESHSET_SET);
  mbw_->add_entities(boundary, this->boundary_vertices);

  return {
    // MOAB's default meshset 0 matches everything.
    {"everywhere", 0},
    {"boundary", boundary}
  };
}
// =============================================================================
std::vector<moab::EntityHandle>
mesh::
compute_boundary_skin_() const
{
  // Find the dimension we're operating on.
  const auto dim =
    this->mbw_->get_number_entities_by_type(0, moab::MBTET) > 0 ? 3 : 2;

  // get all the cell elements on each task
  moab::Range body_entities = this->mbw_->get_entities_by_dimension(0, dim);

  // get face skin
  moab::Skinner tool(this->mbw_->mb.get());
  moab::Range skin_entities;
  moab::ErrorCode rval;
  rval = tool.find_skin(0, body_entities, false, skin_entities);
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);

  // filter out skin elements that are shared with other tasks; they will not
  // be on the true skin
  moab::Range shared_skin;
  rval = this->mcomm_->filter_pstatus(
      skin_entities,
      PSTATUS_SHARED,
      PSTATUS_AND,
      -1,
      &shared_skin
      );
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);

  if (!shared_skin.empty()) {
    skin_entities = subtract(skin_entities, shared_skin);
  }

  // convert range to vector
  std::vector<moab::EntityHandle> skin_entities_vector(
      skin_entities.begin(),
      skin_entities.end()
      );

  return skin_entities_vector;
}
// =============================================================================
moab::Range
mesh::
compute_boundary_vertices_(
    const std::vector<moab::EntityHandle> & boundary_skin
    ) const
{
  // get all vertices on the boundary edges
  const auto verts = this->mbw_->get_adjacencies(
      boundary_skin,
      0,
      false,
      moab::Interface::UNION
      );

  // convert range to vector
  // std::vector<moab::EntityHandle> boundary_verts(verts.begin(), verts.end());

  return verts;
}
// =============================================================================
void
mesh::
mark_subdomains(const std::set<std::shared_ptr<nosh::subdomain>> & subdomains)
{
  const auto & owned_vertices = this->get_owned_vertices();

  for (const auto sd: subdomains) {
    // create meshset
    this->meshsets_[sd->id] = this->mbw_->create_meshset(moab::MESHSET_SET);

    // take care of the vertices
    moab::Range verts = sd->is_boundary_only ?
      this->boundary_vertices :
      owned_vertices;
    for (const auto & vert: verts) {
      const auto x = this->get_coords(vert);
      if (sd->is_inside(x)) {
        mbw_->add_entities(this->meshsets_.at(sd->id), {vert});
      }
    }

    // Take care of the edges.
    // We never need edges on the boundaries, so skip that here.
    if (sd->is_boundary_only) {
      continue;
    }

    const auto edges = this->mbw_->get_entities_by_type(0, moab::MBEDGE);
    for (const auto & edge: edges) {
      // Check if edge midpoint is_inside.
      const auto v = this->mbw_->get_connectivity(edge);
      const auto x0 = this->get_coords(v[0]);
      const auto x1 = this->get_coords(v[1]);
      const auto mp = 0.5 * (x0 + x1);
      if (sd->is_inside(mp)) {
        mbw_->add_entities(this->meshsets_.at(sd->id), {edge});
      }
    }
  }
}
// =============================================================================
moab::Range
mesh::
get_vertices(const std::string & subdomain_id) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      this->meshsets_.count(subdomain_id) == 0,
      "Subdomain \"" << subdomain_id << "\" not found on mesh. "
      << " Did you call mark_subdomains({...}) on the mesh?"
      );

  return this->mbw_->get_entities_by_type(
      this->meshsets_.at(subdomain_id),
      moab::MBVERTEX
      );
}
// =============================================================================
moab::Range
mesh::
get_edges(const std::string & subdomain_id) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
    this->meshsets_.count(subdomain_id) == 0,
    "Subdomain \"" << subdomain_id << "\" not found on mesh. "
      << "Did you call mark_subdomains({...}) on the mesh?"
      );

  return this->mbw_->get_entities_by_type(
      this->meshsets_.at(subdomain_id),
      moab::MBEDGE
      );
}
// =============================================================================
moab::Range
mesh::
get_vertex_tuple(const moab::EntityHandle & edge) const
{
  return this->mbw_->get_adjacencies({edge}, 0, false);
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

  this->mbw_->write_file(
      filename,
      "",
      "PARALLEL=WRITE_PART;"
      );
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
  const moab::Tag tag = this->mbw_->tag_get_handle(tag_name);

  TEUCHOS_ASSERT_EQUALITY(
      this->mbw_->tag_get_data_type(tag),
      moab::DataType::MB_TYPE_DOUBLE
      );

  return this->mbw_->tag_get_data(tag, range);
}
// =============================================================================
Eigen::Vector3d
mesh::
get_coords(
  const moab::EntityHandle vertex
  ) const
{
  return Eigen::Vector3d(this->mbw_->get_coords({vertex}).data());
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh::
get_vector(const std::string & tag_name) const
{
  // get data for all vertices
  moab::Range all_verts = this->mbw_->get_entities_by_dimension(0, 0);

  moab::ErrorCode ierr;
  moab::Range verts;
  ierr = this->mcomm_->filter_pstatus(
      all_verts,
      PSTATUS_NOT_OWNED, PSTATUS_NOT,
      -1,
      &verts
      );
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  auto data = this->get_data(tag_name, verts);

  TEUCHOS_ASSERT_EQUALITY(
    data.size(),
    this->vertices_map_->getNodeNumElements()
    );

  // Set vector values from an existing array (copy)
  return std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(this->vertices_map_),
      Teuchos::ArrayView<double>(data)
      );
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh::
get_complex_vector(const std::string & tag_name) const
{
  // get data for all vertices
  moab::Range all_verts = this->mbw_->get_entities_by_dimension(0, 0);

  moab::ErrorCode ierr;
  moab::Range verts;
  ierr = this->mcomm_->filter_pstatus(
      all_verts,
      PSTATUS_NOT_OWNED, PSTATUS_NOT,
      -1,
      &verts
      );
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  auto data = this->get_data(tag_name, verts);

  TEUCHOS_ASSERT_EQUALITY(
    data.size(),
    this->complex_map_->getNodeNumElements()
    );

  // Set vector values from an existing array (copy)
  return std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(this->complex_map_),
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
  const moab::Range verts = this->mbw_->get_entities_by_dimension(0, 0);

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
    const std::string & name
    ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*write_time_);
#endif

  // get/create handle
  const auto out = this->mbw_->tag_get_handle(
      name,
      1,
      moab::MB_TYPE_DOUBLE,
      moab::MB_TAG_EXCL | moab::MB_TAG_DENSE
      );

  const moab::Tag handle = std::get<0>(out);
  //const bool created = std::get<1>(out);

  // get vertices for which to set the data
  moab::Range all_verts = this->mbw_->get_entities_by_dimension(0, 0);
  moab::ErrorCode ierr;
  moab::Range verts;
  ierr = this->mcomm_->filter_pstatus(
      all_verts,
      PSTATUS_NOT_OWNED, PSTATUS_NOT,
      -1,
      &verts
      );
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  // get data
  const auto data = x.getData();

  TEUCHOS_ASSERT_EQUALITY(
    data.size(),
    this->vertices_map_->getNodeNumElements()
    );

  // set data
  mbw_->tag_set_data(
      handle,
      verts,
      data.get()
      );

  return;
}
// =============================================================================
const std::vector<Teuchos::Tuple<int,2>>
mesh::
build_edge_lids_() const
{
  const std::vector<edge> edges = this->my_edges();

  std::vector<Teuchos::Tuple<int,2>> _edge_lids(edges.size());

  for (std::size_t k = 0; k < edges.size(); k++) {
    _edge_lids[k] = Teuchos::tuple(
        (int)this->local_index(std::get<0>(edges[k])),
        (int)this->local_index(std::get<1>(edges[k]))
        );
  }

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
const std::vector<Teuchos::Tuple<int,2>>
mesh::
build_edge_gids_() const
{
  const std::vector<edge> edges = this->my_edges();

  std::vector<Teuchos::Tuple<int,2>> _edge_gids(edges.size());

  const moab::Tag gid = this->mbw_->tag_get_handle("GLOBAL_ID");

  for (std::size_t k = 0; k < edges.size(); k++) {
    // get the global IDs of the vertices
    const std::vector<moab::EntityHandle> entities =
      {std::get<0>(edges[k]), std::get<1>(edges[k])};
    const auto global_ids = this->mbw_->tag_get_int_data(gid, entities);

    _edge_gids[k] = Teuchos::tuple(
        global_ids[0],
        global_ids[1]
        );
  }

  return _edge_gids;
}
// =============================================================================
const std::vector<Teuchos::Tuple<int,4>>
mesh::
build_edge_gids_complex_() const
{
  const std::vector<edge> edges = this->my_edges();

  std::vector<Teuchos::Tuple<int,4>> _edge_gids_complex(edges.size());

  const moab::Tag gid = this->mbw_->tag_get_handle("GLOBAL_ID");

  for (std::size_t k = 0; k < edges.size(); k++) {
    // get the global IDs of the vertices
    const std::vector<moab::EntityHandle> entities =
      {std::get<0>(edges[k]), std::get<1>(edges[k])};
    const auto global_ids = this->mbw_->tag_get_int_data(gid, entities);

    _edge_gids_complex[k] =
      Teuchos::tuple(
          2*global_ids[0], 2*global_ids[0] + 1,
          2*global_ids[1], 2*global_ids[1] + 1
          );
  }

  return _edge_gids_complex;
}
// =============================================================================
moab::Range
mesh::
get_owned_vertices() const
{
  // TODO resurrect?
  //const auto mb = this->mcomm_->get_moab();

  // get all entities
  moab::Range all_verts = this->mbw_->get_entities_by_dimension(0, 0);

  // filter out only owned
  moab::ErrorCode ierr;
  moab::Range verts;
  ierr = this->mcomm_->filter_pstatus(
      all_verts,
      PSTATUS_NOT_OWNED, PSTATUS_NOT,
      -1,
      &verts
      );
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  return verts;
}
// =============================================================================
const std::vector<int>
mesh::
get_owned_gids_() const
{
  moab::Tag gid = mbw_->tag_get_handle("GLOBAL_ID");
  return this->mbw_->tag_get_int_data(gid, this->get_owned_vertices());
}
// =============================================================================
const std::vector<int>
mesh::
get_overlap_gids_() const
{
  // TODO remove?
  //const auto mb = this->mcomm_->get_moab();

  // get owned
  moab::Range all = this->mbw_->get_entities_by_dimension(0, 0);

  // Get entities shared with all other processors
  moab::Range shared;
  moab::ErrorCode ierr;
  ierr = mcomm_->get_shared_entities(-1, shared, 0); // only verts
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  // merge
  all.merge(shared);

  // get the corresponding global IDs
  const moab::Tag gid = mbw_->tag_get_handle("GLOBAL_ID");

  return this->mbw_->tag_get_int_data(gid, all);
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
  // Given all IDs, the base is actually redundant. It has to be <= the
  // minimal global index, and is often taken to be 0. However, there are some
  // bugs that manifest if the minimal global index does not equal the base,
  // e.g., <https://github.com/trilinos/Trilinos/issues/70>.
  // Hence, set it to 1 as dictated by MOAB.
  //
  // The will fail for complexified maps where the minimal ID is 2.
  // TODO derive the base from all ID vectors.
  const int base = 1;
  return std::make_shared<Tpetra::Map<int,int>>(
      -1,
      ids,
      base,
      Teuchos::rcp(this->comm)
      );
}
// =============================================================================
mesh::entity_relations
mesh::
build_entity_relations_()
{
  // get the number of 3D entities
  const int num3d = this->mbw_->get_number_entities_by_dimension(0, 3);

  const int dim = (num3d > 0) ? 3 : 2;

  // Get regions, by dimension, so we stay generic to entity type
  const moab::Range elems = this->mbw_->get_entities_by_dimension(0, dim);

  // get and create all edges adjacent to cells
  const moab::Range edges = this->mbw_->get_adjacencies(
      elems,
      1,
      true,
      moab::Interface::UNION
      );

  // create cell->edge relation
  std::vector<std::vector<moab::EntityHandle>> cell_edges(elems.size());
  for (size_t k = 0; k < elems.size(); k++) {
    const auto ce = this->mbw_->get_adjacencies(
          {elems[k]},
          1,
          true,
          moab::Interface::UNION
          );

    cell_edges[k].resize(ce.size());
    for (size_t i = 0; i < ce.size(); i++) {
      cell_edges[k][i] = ce[i];
    }
  }

  // create edge->node relation
  std::vector<std::tuple<moab::EntityHandle, moab::EntityHandle>>
    edge_vertices(edges.size());

  for (size_t k = 0; k < edges.size(); k++) {
    moab::Range verts = this->mbw_->get_adjacencies(
        {edges[k]},
        0,
        true,
        moab::Interface::UNION
        );

    edge_vertices[k] = std::make_tuple(verts[0], verts[1]);
  }

  mesh::entity_relations relations = {edge_vertices, cell_edges};

  //// for testing
  //moab::Range verts = this->mbw_->mb->get_adjacencies(
  //    elems,
  //    0,
  //    true,
  //    moab::Interface::UNION
  //    );

  return relations;
}
// =============================================================================
Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
mesh::
build_graph() const
{
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from
  // the owned_vertices map, and the overlapping one from the
  // overlap_vertices.
  // Let's illustrate the implications with the example of the matrix
  //   [ 2 1   ]
  //   [ 1 2 1 ]
  //   [   1 2 ].
  // Suppose subdomain 1 consists of node 1, subdomain 2 of node 3,
  // and node 2 forms the boundary between them.
  // For two processes, if process 1 owns vertices 1 and 2, the matrix
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
  // In the general case, assuming that the number of vertices adjacent
  // to a boundary (on one side) are approximately the number of
  // vertices on that boundary, there is not much difference in
  // communication between the patterns.
  // What does differ, though, is the workload on the processes
  // during the computation phase: Process 1 that owns the whole
  // boundary, has to compute more than process 2.
  // Notice, however, that the total number of computations is
  // lower in scenario 1 (7 vs. 8 FLOPs); the same is true for
  // storage.
  // Hence, it comes down to the question whether or not the
  // mesh generator provided a fair share of the boundary vertices.
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
  //auto graph = Teuchos::rcp(new Tpetra::CrsGraph<int, int>(
  //    Teuchos::rcp(nonoverlap_map),
  //    Teuchos::rcp(nonoverlap_map),
  //    0
  //    ));
  const auto graph = Tpetra::createCrsGraph(Teuchos::rcp(nonoverlap_map));

  const std::vector<edge> edges = this->my_edges();

  // Loop over all edges and put entries wherever two vertices are connected.
  // TODO check if we can use LIDs here
  for (size_t k = 0; k < edges.size(); k++) {
    const Teuchos::Tuple<int,2> & idx = this->edge_gids[k];
    for (int i = 0; i < 2; i++) {
      graph->insertGlobalIndices(idx[i], idx);
      // graph->insertLocalIndices(idx[i], idx);
    }
  }

  graph->fillComplete();

  return graph;
}
// =============================================================================
Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
mesh::
build_complex_graph() const
{
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from
  // the owned_vertices map, and the overlapping one from the
  // overlap_vertices.
  // Let's illustrate the implications with the example of the matrix
  //   [ 2 1   ]
  //   [ 1 2 1 ]
  //   [   1 2 ].
  // Suppose subdomain 1 consists of node 1, subdomain 2 of node 3,
  // and node 2 forms the boundary between them.
  // For two processes, if process 1 owns vertices 1 and 2, the matrix
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
  // In the general case, assuming that the number of vertices adjacent
  // to a boundary (on one side) are approximately the number of
  // vertices on that boundary, there is not much difference in
  // communication between the patterns.
  // What does differ, though, is the workload on the processes
  // during the computation phase: Process 1 that owns the whole
  // boundary, has to compute more than process 2.
  // Notice, however, that the total number of computations is
  // lower in scenario 1 (7 vs. 8 FLOPs); the same is true for
  // storage.
  // Hence, it comes down to the question whether or not the
  // mesh generator provided a fair share of the boundary vertices.
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
  // const auto graph = Teuchos::rcp(new Tpetra::CrsGraph<int, int>(
  //     Teuchos::rcp(nonoverlap_map),
  //     Teuchos::rcp(nonoverlap_map),
  //     0
  //     ));
  const auto graph = Tpetra::createCrsGraph(Teuchos::rcp(nonoverlap_map));

  const std::vector<edge> edges = this->my_edges();

  // Loop over all edges and put entries wherever two vertices are connected.
  // TODO check if we can use LIDs here
  for (size_t k = 0; k < edges.size(); k++) {
    const Teuchos::Tuple<int,4> & idx = this->edge_gids_complex[k];
    for (int i = 0; i < 4; i++) {
#ifdef NDEBUG
      TEUCHOS_ASSERT_INEQUALITY(idx[i], >=, 0);
#endif
      graph->insertGlobalIndices(idx[i], idx);
      //graph->insertLocalIndices(idx[i], idx);
    }
  }

  graph->fillComplete();

  return graph;
}
// =============================================================================
Eigen::Vector3d
mesh::
compute_triangle_circumcenter_(
  const std::vector<Eigen::Vector3d> &vertices
  ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(vertices.size(), 3);
#endif
  return this->compute_triangle_circumcenter_(vertices[0], vertices[1], vertices[2]);
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
      "It seems that the vertices \n"
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
std::vector<double>
mesh::
compute_triangle_splitting_(
    const std::vector<moab::EntityHandle> & conn
    ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(conn.size(), 3);
#endif

  std::vector<double> splitting = {0.0, 0.0, 0.0};

  // Fetch the nodal positions into 'local_node_coords'.
  std::vector<double> coords = this->mbw_->get_coords(conn);

  std::vector<Eigen::Vector3d> local_node_coords(conn.size());
  for (size_t i = 0; i < conn.size(); i++) {
    // TODO do something smarter than copying here
    local_node_coords[i][0] = coords[3*i];
    local_node_coords[i][1] = coords[3*i + 1];
    local_node_coords[i][2] = coords[3*i + 2];
  }

  // compute the circumcenter of the cell
  const Eigen::Vector3d cc =
    compute_triangle_circumcenter_(local_node_coords);

  // Iterate over the edges (aka pairs of vertices).
  for (size_t e0 = 0; e0 < 3; e0++) {
    const Eigen::Vector3d &x0 = local_node_coords[e0];
    for (size_t e1 = e0+1; e1 < 3; e1++) {
      const Eigen::Vector3d &x1 = local_node_coords[e1];
      // Get the other node.
      const unsigned int other = this->get_other_index_(e0, e1);

      double edge_length = (x1-x0).norm();

      // Compute the (n-1)-dimensional covolume.
      const Eigen::Vector3d &other0 = local_node_coords[other];
      double covolume = this->compute_covolume2d_(cc, x0, x1, other0);
      // Compute the contributions to the finite volumes of the adjacent
      // edges.
      double pyramid_volume = 0.5 * edge_length * covolume / 2;
      splitting[e0] += pyramid_volume;
      splitting[e1] += pyramid_volume;
    }
  }

  return splitting;
}
// =============================================================================
unsigned int
mesh::
get_other_index_(unsigned int e0, unsigned int e1) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_INEQUALITY(e0, !=, e1);
#endif

  // Get the index in [0,1,2] which is not e0, e1.
  if (0 != e0 && 0 != e1)
    return 0;
  else if (1 != e0 && 1 != e1)
    return 1;
  else if (2 != e0 && 2 != e1)
    return 2;
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        "illegal"
        );
}
// =============================================================================
double
mesh::
compute_covolume2d_(
    const Eigen::Vector3d &cc,
    const Eigen::Vector3d &x0,
    const Eigen::Vector3d &x1,
    const Eigen::Vector3d &other0
    ) const
{
  // edge midpoint
  Eigen::Vector3d mp = 0.5 * (x0 + x1);

  double coedge_length = (mp - cc).norm();

  // The only difficulty here is to determine whether the length of coedge is
  // to be taken positive or negative.
  // To this end, make sure that the order (x0, cc, mp) is of the same
  // orientation as (x0, other0, mp).
  Eigen::Vector3d cell_normal = (other0 - x0).cross(mp - x0);
  Eigen::Vector3d cc_normal = (cc - x0).cross(mp - x0);

  // copysign takes the absolute value of the first argument and the sign of
  // the second.
  return copysign(coedge_length, cc_normal.dot(cell_normal));
}
// =============================================================================
}  // namespace nosh
