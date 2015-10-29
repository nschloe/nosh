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
#include "mesh_tri.hpp"

#include <vector>
#include <set>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <stk_mesh/base/MetaData.hpp>

namespace nosh
{
// =============================================================================
mesh_tri::
mesh_tri(
    const std::shared_ptr<const Teuchos::Comm<int>> & _comm,
    const std::shared_ptr<stk::io::StkMeshIoBroker> & broker,
      const std::set<std::string> allocated_vector_names
    ) :
  mesh(_comm, broker, allocated_vector_names),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  compute_edge_coefficients_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tri::compute_edge_coefficients"
        )),
  compute_control_volumes_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tri::compute_control_volumes"
        )),
  compute_boundary_nodes_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tri::compute_boundary_nodes"
        )),
#endif
  control_volumes_(this->compute_control_volumes_()),
  edge_coefficients_(this->compute_edge_coefficients_()),
  boundary_nodes_(this->compute_boundary_nodes_())
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

  std::vector<stk::mesh::Entity> cells = this->get_owned_cells();
  unsigned int num_cells = cells.size();

  size_t num_eges = edge_data_.edge_nodes.size();

  std::vector<double> _edge_coefficients(num_eges);

  const vector_fieldType & coords_field = get_node_field("coordinates");

  // Calculate the contributions edge by edge.
  for (unsigned int k = 0; k < num_cells; k++) {
    // Get edge coordinates.
    size_t num_local_edges = edge_data_.cell_edges[k].size();
    std::vector<Eigen::Vector3d> local_edge_coords(num_local_edges);
    for (size_t i = 0; i < num_local_edges; i++) {
      const edge & e = edge_data_.edge_nodes[edge_data_.cell_edges[k][i]];
      local_edge_coords[i] =
        this->get_node_value(coords_field, std::get<1>(e))
        - this->get_node_value(coords_field, std::get<0>(e));
    }

    Eigen::VectorXd edge_coeffs =
      edge_coefficients_numerically_(local_edge_coords);

    // Fill the edge coefficients into the vector.
    for (size_t i = 0; i < num_local_edges; i++) {
      _edge_coefficients[edge_data_.cell_edges[k][i]] += edge_coeffs[i];
    }
  }

  return _edge_coefficients;
}
// =============================================================================
Eigen::VectorXd
mesh_tri::
edge_coefficients_numerically_(
  const std::vector<Eigen::Vector3d> edges
  ) const
{
  size_t num_eges = edges.size();

  // Build an equation system for the edge coefficients alpha_k.
  // They fulfill
  //
  //    |simplex| * <u,v> = \sum_{edges e_i} alpha_i <u,e_i> <e_i,v>
  //
  // for any pair of vectors u, v in the plane of the triangle.
  //
  const double vol = 0.5 * (edges[0].cross(edges[1])).norm();

  Eigen::MatrixXd A(num_eges, num_eges);
  Eigen::VectorXd rhs(num_eges);

  // Build the equation system:
  // The equation
  //
  //    |simplex| ||u||^2 = \sum_i \alpha_i <u,e_i> <e_i,u>
  //
  // has to hold for all vectors u in the plane spanned by the edges,
  // particularly by the edges themselves.
  //
  // Only fill the upper part of the Hermitian matrix.
  //
  for (size_t i = 0; i < num_eges; i++) {
    double alpha = edges[i].dot(edges[i]);
    rhs(i) = vol * alpha;
    A(i,i) = alpha * alpha;
    for (size_t j = i+1; j < num_eges; j++) {
      A(i, j) = edges[i].dot(edges[j]) * edges[j].dot(edges[i]);
      A(j, i) = A(i, j);
    }
  }

  // Solve the equation system for the alpha_i.  The system is symmetric and,
  // if the simplex is not degenerate, positive definite.
  //return A.ldlt().solve(rhs);
  return A.fullPivLu().solve(rhs);
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
  TEUCHOS_ASSERT(io_broker_);
#endif

  // Create temporaries to hold the overlap values for control volumes.
  Tpetra::Vector<double,int,int> cv_overlap(Teuchos::rcp(nodes_overlap_map_));

  this->compute_control_volumes_t_(cv_overlap);

  // Export control volumes to a non-overlapping map, and sum the entries.
  Teuchos::RCP<const Tpetra::Export<int,int>> exporter = Tpetra::createExport(
      Teuchos::rcp(nodes_overlap_map_),
      Teuchos::rcp(nodes_map_)
      );

  auto _control_volumes = std::make_shared<Tpetra::Vector<double,int,int>>(
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
  std::vector<stk::mesh::Entity> cells = this->get_owned_cells();
  size_t num_cells = cells.size();

  Teuchos::ArrayRCP<double> cv_data = cv_overlap.getDataNonConst();

  const vector_fieldType & coords_field = get_node_field("coordinates");

  // Calculate the contributions to the finite volumes cell by cell.
  for (size_t k = 0; k < num_cells; k++) {
    const stk::mesh::Entity * local_nodes =
      io_broker_->bulk_data().begin_nodes(cells[k]);
    unsigned int num_local_nodes = io_broker_->bulk_data().num_nodes(cells[k]);

#ifndef NDEBUG
    // Confirm that we always have the same simplices.
    TEUCHOS_ASSERT_EQUALITY(num_local_nodes, 3);
#endif

    // Fetch the nodal positions into 'local_nodes'.
    //const std::vector<Eigen::Vector3d> local_node_coords =
    //this->get_nodeCoordinates_(local_nodes);
    std::vector<Eigen::Vector3d> local_node_coords(num_local_nodes);
    for (unsigned int i = 0; i < num_local_nodes; i++) {
      local_node_coords[i] = this->get_node_value(coords_field, local_nodes[i]);
    }

    // compute the circumcenter of the cell
    const Eigen::Vector3d cc = compute_triangle_circumcenter_(local_node_coords);

    // Iterate over the edges.
    // As true edge entities are not available here, loop over all pairs of
    // local nodes.
    for (unsigned int e0 = 0; e0 < num_local_nodes; e0++) {
      const Eigen::Vector3d &x0 = local_node_coords[e0];
      const int gid0 = this->gid(local_nodes[e0]);
      const int lid0 = nodes_overlap_map_->getLocalElement(gid0);
#ifndef NDEBUG
      TEUCHOS_ASSERT_INEQUALITY(lid0, >=, 0);
#endif
      for (unsigned int e1 = e0+1; e1 < num_local_nodes; e1++) {
        const Eigen::Vector3d &x1 = local_node_coords[e1];
        const int gid1 = this->gid(local_nodes[e1]);
        const int lid1 = nodes_overlap_map_->getLocalElement(gid1);
#ifndef NDEBUG
        TEUCHOS_ASSERT_INEQUALITY(lid1, >=, 0);
#endif
        // Get the other node.
        const unsigned int other = this->get_other_index_(e0, e1);

        double edge_length = (x1-x0).norm();

        // Compute the (n-1)-dimensional covolume.
        double covolume;
        const Eigen::Vector3d &other0 = local_node_coords[other];
        covolume = this->compute_covolume2d_(cc, x0, x1, other0);
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
        cv_data[lid0] += pyramid_volume;
        cv_data[lid1] += pyramid_volume;
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
std::set<stk::mesh::Entity>
mesh_tri::
compute_boundary_nodes_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*compute_boundary_nodes_time_);
#endif

  auto _my_edges = this->get_overlap_edges();

  std::set<stk::mesh::Entity> _boundary_nodes;
  for (size_t k = 0; k < _my_edges.size(); k++) {
    // if the edge has one element, it's on the boundary
    if (io_broker_->bulk_data().num_elements(_my_edges[k]) == 1) {
      stk::mesh::Entity const * nodes =
        io_broker_->bulk_data().begin_nodes(_my_edges[k]);
#ifndef NDEBUG
      TEUCHOS_ASSERT_EQUALITY(io_broker_->bulk_data().num_nodes(_my_edges[k]), 2);
#endif
      _boundary_nodes.insert(nodes[0]);
      _boundary_nodes.insert(nodes[1]);
    }
  }

  return _boundary_nodes;
}
// =============================================================================
}  // namespace nosh
