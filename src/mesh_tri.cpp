#include "mesh_tri.hpp"

#include <vector>
#include <set>

#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <moab/Skinner.hpp>
#include <MBParallelConventions.h>

namespace nosh
{
// =============================================================================
mesh_tri::
mesh_tri(
    const std::shared_ptr<const Teuchos::Comm<int>> & _comm,
    const std::shared_ptr<moab::ParallelComm> & mcomm,
    const std::shared_ptr<moab::Core> & mb
    ) :
  mesh(_comm, mcomm, mb)
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  ,compute_edge_coefficients_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tri::compute_edge_coefficients"
        ))
  ,compute_control_volumes_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tri::compute_control_volumes"
        ))
  ,compute_boundary_nodes_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tri::compute_boundary_nodes"
        ))
#endif
  ,control_volumes_(this->compute_control_volumes_())
  ,edge_coefficients_(this->compute_edge_coefficients_())
  ,boundary_nodes_(this->compute_boundary_nodes_())
{
}
// =============================================================================
mesh_tri::
~mesh_tri()
{
}
// =============================================================================
std::vector<double>
mesh_tri::
compute_edge_coefficients_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*compute_edge_coefficients_time_);
#endif

  moab::Range cells = this->mbw_->get_entities_by_dimension(0, 2);

  size_t num_cells = cells.size();

  size_t num_edges = relations_.edge_nodes.size();

  // compute all coordinates
  std::vector<Eigen::Vector3d> edge_coords(num_edges);
  for (size_t k = 0; k < num_edges; k++) {
    auto tmp1 = std::get<0>(relations_.edge_nodes[k]);
    std::vector<double> coords0 = this->mbw_->get_coords({tmp1});

    tmp1 = std::get<1>(relations_.edge_nodes[k]);
    std::vector<double> coords1 = this->mbw_->get_coords({tmp1});

    edge_coords[k][0] = coords0[0] - coords1[0];
    edge_coords[k][1] = coords0[1] - coords1[1];
    edge_coords[k][2] = coords0[2] - coords1[2];
  }

  std::vector<double> _edge_coefficients(num_edges);

  // Compute the contributions cell by cell.
  for (size_t k = 0; k < num_cells; k++) {
    const std::vector<size_t> edge_idxs = {
      this->local_index(relations_.cell_edges[k][0]),
      this->local_index(relations_.cell_edges[k][1]),
      this->local_index(relations_.cell_edges[k][2])
    };
    const std::vector<Eigen::Vector3d> local_edge_coords = {
      edge_coords[edge_idxs[0]],
      edge_coords[edge_idxs[1]],
      edge_coords[edge_idxs[2]]
    };
    // const moab::EntityHandle * conn = NULL;
    // int numV = 0;
    // ierr = this->mbw_->mb->get_connectivity(cells[k], conn, numV);
    // TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

    // std::vector<double> coords(3 * numV);
    // ierr = this->mbw_->mb->get_coords(conn, numV, &coords[0]);
    // TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

    // // Get edge coordinates.
    // const unsigned int num_local_edges = (numV - 1) * numV / 2;
    // std::vector<Eigen::Vector3d> local_edge_coords(num_local_edges);
    // int k0 = 0;
    // for (int i = 0; i < numV; i++) {
    //   for (int j = i+1; j < numV; j++) {
    //     local_edge_coords[k0][0] = coords[3*i] - coords[3*j];
    //     local_edge_coords[k0][1] = coords[3*i + 1] - coords[3*j + 1];
    //     local_edge_coords[k0][2] = coords[3*i + 2] - coords[3*j + 2];
    //     k0++;
    //   }
    // }

    Eigen::VectorXd edge_coeffs =
      edge_coefficients_numerically_(local_edge_coords);

    // Fill the edge coefficients into the vector.
    for (int i = 0; i < edge_coeffs.size(); i++) {
      const size_t edge_idx = this->local_index(relations_.cell_edges[k][i]);
      // const int edge_id = this->liddd
      _edge_coefficients[edge_idx] += edge_coeffs[i];
    }
  }

  return _edge_coefficients;
}
// =============================================================================
Eigen::VectorXd
mesh_tri::
edge_coefficients_numerically_(
  const std::vector<Eigen::Vector3d> & edges
  ) const
{
  size_t num_edges = edges.size();
  TEUCHOS_ASSERT_EQUALITY(num_edges, 3);

  // Build an equation system for the edge coefficients alpha_k.
  // They fulfill
  //
  //    |simplex| * <u,v> = \sum_{edges e_i} alpha_i <u,e_i> <e_i,v>
  //
  // for any pair of vectors u, v in the plane of the triangle.
  //
  const double vol = 0.5 * (edges[0].cross(edges[1])).norm();

  Eigen::Matrix3d A;
  Eigen::Vector3d rhs;

  // Build the equation system:
  // The equation
  //
  //    |simplex| ||u||^2 = \sum_i \alpha_i <u,e_i> <e_i,u>
  //
  // has to hold for all vectors u in the plane spanned by the edges,
  // particularly by the edges themselves.
  //
  for (size_t i = 0; i < num_edges; i++) {
    double alpha = edges[i].dot(edges[i]);
    rhs(i) = vol * alpha;
    A(i, i) = alpha * alpha;
    for (size_t j = i+1; j < num_edges; j++) {
      A(i, j) = edges[i].dot(edges[j]) * edges[j].dot(edges[i]);
      A(j, i) = A(i, j);
    }
  }

  // Solve the equation system for the alpha_i.  The system is symmetric and,
  // if the simplex is not degenerate, positive definite.
  //return A.ldlt().solve(rhs);
  const auto x = A.fullPivLu().solve(rhs);
  //auto x = A.colPivHouseholderQr().solve(rhs);

  return x;
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh_tri::
compute_control_volumes_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*compute_control_volumes_time_);
#endif
#ifndef NDEBUG
  TEUCHOS_ASSERT(nodes_map_);
  TEUCHOS_ASSERT(nodes_overlap_map_);
#endif

  auto _control_volumes = std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(nodes_map_)
      );

  // Create temporaries to hold the overlap values for control volumes.
  Tpetra::Vector<double,int,int> cv_overlap(Teuchos::rcp(nodes_overlap_map_));

  this->compute_control_volumes_t_(cv_overlap);

  // Export control volumes to a non-overlapping map, and sum the entries.
  Teuchos::RCP<const Tpetra::Export<int,int>> exporter = Tpetra::createExport(
      Teuchos::rcp(nodes_overlap_map_),
      Teuchos::rcp(nodes_map_)
      );

  _control_volumes->doExport(cv_overlap, *exporter, Tpetra::ADD);

  return _control_volumes;
}
// =============================================================================
void
mesh_tri::
compute_control_volumes_t_(Tpetra::Vector<double,int,int> & cv_overlap) const
{
  // get owned entities
  moab::Range cells = this->mbw_->get_entities_by_dimension(0, 2);

  Teuchos::ArrayRCP<double> cv_data = cv_overlap.getDataNonConst();
  //const vector_fieldType & coords_field = get_node_field("coordinates");

  // std::vector<double> coords;
  // ierr = mcomm_->getMoab()->get_vertex_coordinates(coords);
  // TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  // Calculate the contributions to the finite volumes cell by cell.
  for (size_t k = 0; k < cells.size(); k++) {
    const auto conn = this->mbw_->get_connectivity(cells[k]);

#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(conn.size(), 3);
#endif

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

    // Iterate over the edges (aka pairs of nodes).
    for (size_t e0 = 0; e0 < conn.size(); e0++) {
      const Eigen::Vector3d &x0 = local_node_coords[e0];
      for (size_t e1 = e0+1; e1 < conn.size(); e1++) {
        const Eigen::Vector3d &x1 = local_node_coords[e1];
        // Get the other node.
        const unsigned int other = this->get_other_index_(e0, e1);

        double edge_length = (x1-x0).norm();

        // Compute the (n-1)-dimensional covolume.
        const Eigen::Vector3d &other0 = local_node_coords[other];
        double covolume = this->compute_covolume2d_(cc, x0, x1, other0);
        // The problem with counting the average thickness in 2D is the
        // following.  Ideally, one would want to loop over all edges, add the
        // midpoint value of the thickness to both of the edge end points, and
        // eventually loop over all endpoints and divide by the number of edges
        // (connections) they have with neighboring nodes).
        // Unfortunately, this is impossible now b/c there's no edge generation
        // for shells in Trilinos yet (2011-04-15).
        // As a workaround, one could loop over all cells, and then all pairs
        // of nodes to retrieve the edges. In 2D, almost all of the edges would
        // be counted twice this way as they belong to two cells. This is true
        // for all but the boundary edges. Again, it is difficult (impossible?)
        // to know what the boundary edges are, and hence which values to
        // divide by 2. Dividing them all by two would result in an
        // artificially lower thickness near the boundaries. This is not what
        // we want.

        // Compute the contributions to the finite volumes of the adjacent
        // edges.
        double pyramid_volume = 0.5 * edge_length * covolume / 2;
        // The EntityHandle (conn) is a local identifier, MOAB indices are
        // 1-based.
        cv_data[conn[e0] - 1] += pyramid_volume;
        cv_data[conn[e1] - 1] += pyramid_volume;
      }
    }
  }

  return;
}
// =============================================================================
double
mesh_tri::
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
unsigned int
mesh_tri::
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
std::set<moab::EntityHandle>
mesh_tri::
compute_boundary_nodes_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*compute_boundary_nodes_time_);
#endif

  // get all the cell elements on each task
  moab::Range cells = this->mbw_->get_entities_by_dimension(0, 2);

  // get face skin
  moab::Skinner tool(this->mbw_->mb.get());
  moab::Range edges;
  moab::ErrorCode rval;
  rval = tool.find_skin(0, cells, false, edges);
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);

  // filter out edges that are shared with other tasks; they will not be on the
  // true skin
  moab::Range shared_edges;
  rval = this->mcomm_->filter_pstatus(
      edges,
      PSTATUS_SHARED,
      PSTATUS_AND,
      -1,
      &shared_edges
      );
  TEUCHOS_ASSERT_EQUALITY(rval, moab::MB_SUCCESS);

  if (!shared_edges.empty()) {
    edges = subtract(edges, shared_edges);
  }

  // get all vertices on the remaining edges
  const auto verts = this->mbw_->get_adjacencies(
      edges,
      0,
      false,
      moab::Interface::UNION
      );

  // convert range to set
  // TODO perhaps there is better way?
  std::set<moab::EntityHandle> boundary_verts;
  for (size_t k = 0; k < verts.size(); k++) {
    boundary_verts.insert(verts[k]);
  }

  return boundary_verts;
}
// =============================================================================
}  // namespace nosh
