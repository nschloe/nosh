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
#include "mesh_tetra.hpp"

#include <memory>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include <Eigen/Dense>

namespace nosh
{
// =============================================================================
mesh_tetra::
mesh_tetra(
    const std::shared_ptr<const Teuchos::Comm<int>> & _comm,
    const std::shared_ptr<moab::ParallelComm> & mcomm,
    const std::shared_ptr<moab::Core> & mb
    ) :
  mesh(_comm, mcomm, mb),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  compute_edge_coefficients_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tetra::compute_edge_coefficients"
        )),
  compute_control_volumes_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tetra::compute_control_volumes"
        )),
  compute_boundary_nodes_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tetra::compute_boundary_nodes"
        ))
#endif
  ,control_volumes_(this->compute_control_volumes_())
  ,edge_coefficients_(this->compute_edge_coefficients_())
  ,boundary_nodes_(this->compute_boundary_nodes_())
{
}
// =============================================================================
mesh_tetra::
~mesh_tetra()
{
}
// =============================================================================
std::vector<double>
mesh_tetra::
compute_edge_coefficients_() const
{
  std::cout << ">> compute_edge_coefficients_" << std::endl;
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*compute_edge_coefficients_time_);
#endif
  moab::ErrorCode ierr;

  moab::Range cells;
  ierr = this->mb_->get_entities_by_dimension(0, 3, cells);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  size_t num_cells = cells.size();

  size_t num_edges = edge_data_.edge_nodes.size();

  // compute all coordinates
  std::vector<Eigen::Vector3d> edge_coords(num_edges);
  for (size_t k = 0; k < num_edges; k++) {
    std::vector<double> coords0(3);
    auto tmp1 = std::get<0>(edge_data_.edge_nodes[k]);
    ierr = this->mb_->get_coords(&tmp1, 1, &coords0[0]);
    TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

    std::vector<double> coords1(3);
    tmp1 = std::get<1>(edge_data_.edge_nodes[k]);
    ierr = this->mb_->get_coords(&tmp1, 1, &coords1[0]);
    TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

    edge_coords[k][0] = coords0[0] - coords1[0];
    edge_coords[k][1] = coords0[1] - coords1[1];
    edge_coords[k][2] = coords0[2] - coords1[2];
  }

  std::vector<double> _edge_coefficients(num_edges);

  // Compute the contributions edge by edge.
  for (size_t k = 0; k < num_cells; k++) {
    const std::vector<size_t> edge_idxs = {
      this->local_index(edge_data_.cell_edges[k][0]),
      this->local_index(edge_data_.cell_edges[k][1]),
      this->local_index(edge_data_.cell_edges[k][2]),
      this->local_index(edge_data_.cell_edges[k][3]),
      this->local_index(edge_data_.cell_edges[k][4]),
      this->local_index(edge_data_.cell_edges[k][5])
    };

    const std::vector<Eigen::Vector3d> local_edge_coords = {
      edge_coords[edge_idxs[0]],
      edge_coords[edge_idxs[1]],
      edge_coords[edge_idxs[2]],
      edge_coords[edge_idxs[3]],
      edge_coords[edge_idxs[4]],
      edge_coords[edge_idxs[5]]
    };

    auto edge_coeffs = edge_coefficients_cell_(local_edge_coords);

    // Fill the edge coefficients into the vector.
    for (int i = 0; i < edge_coeffs.size(); i++) {
      const size_t edge_idx = this->local_index(edge_data_.cell_edges[k][i]);
      _edge_coefficients[edge_idx] += edge_coeffs[i];
    }
  }

  std::cout << "   compute_edge_coefficients_ >>" << std::endl;
  return _edge_coefficients;
}
// =============================================================================
Eigen::VectorXd
mesh_tetra::
edge_coefficients_cell_(
  const std::vector<Eigen::Vector3d> & edges
  ) const
{
  size_t num_edges = edges.size();

  // Build an equation system for the edge coefficients alpha_k.
  // They fulfill
  //
  //    |simplex| * <u,v> = \sum_{edges e_i} alpha_i <u,e_i> <e_i,v>
  //
  // for any pair of vectors u, v in the plane of the triangle.
  //
  double vol;
  // TODO Come up with a cleaner solution here.
  try {
    vol = get_tetrahedron_volume_(edges[0], edges[1], edges[2]);
  } catch(...) {
    // If computing the volume throws an exception, then the edges chosen
    // happened to be conplanar. Changing one of those fixes this.
    vol = get_tetrahedron_volume_(edges[0], edges[1], edges[3]);
  }

  Eigen::Matrix<double, 6, 6> A;
  Eigen::VectorXd rhs(6);

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
  const auto x = A.fullPivLu().solve(rhs);
  //const auto x = A.colPivHouseholderQr().solve(rhs);
  //const auto x = A.ldlt().solve(rhs);

  return x;
}
// =============================================================================
std::shared_ptr<Tpetra::Vector<double,int,int>>
mesh_tetra::
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

  // Create temporaries to hold the overlap values for control volumes and
  // average thickness.
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
mesh_tetra::
compute_control_volumes_t_(Tpetra::Vector<double,int,int> & cv_overlap) const
{
  // get owned entities
  moab::Range cells;

  moab::ErrorCode ierr;

  ierr = this->mb_->get_entities_by_dimension(0, 3, cells);
  TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

  unsigned int num_cells = cells.size();

  Teuchos::ArrayRCP<double> cv_data = cv_overlap.getDataNonConst();

  auto rank = cv_overlap.getMap()->getComm()->getRank();

  // Calculate the contributions to the finite volumes cell by cell.
  for (size_t k = 0; k < num_cells; k++) {
    const moab::EntityHandle * conn = NULL;
    int numV = 0;
    ierr = this->mb_->get_connectivity(cells[k], conn, numV);
    TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(numV, 4);
#endif

    // Fetch the nodal positions into 'local_node_coords'.
    std::vector<double> coords(3 * numV);
    ierr = this->mb_->get_coords(conn, numV, &coords[0]);
    TEUCHOS_ASSERT_EQUALITY(ierr, moab::MB_SUCCESS);

    std::vector<Eigen::Vector3d> local_node_coords(numV);
    for (int i = 0; i < numV; i++) {
      // TODO do something smarter than copying here
      local_node_coords[i][0] = coords[3*i];
      local_node_coords[i][1] = coords[3*i + 1];
      local_node_coords[i][2] = coords[3*i + 2];
    }

    // compute the circumcenter of the cell
    const auto cc = this->compute_tetrahedron_circumcenter_(local_node_coords);

    // Iterate over the edges.
    // As true edge entities are not available here, loop over all pairs of
    // local nodes.
    for (int e0 = 0; e0 < numV; e0++) {
      const Eigen::Vector3d &x0 = local_node_coords[e0];
      for (int e1 = e0+1; e1 < numV; e1++) {
        const Eigen::Vector3d &x1 = local_node_coords[e1];

        // Get the other nodes.
        std::set<unsigned int> other_set = this->get_other_indices_(e0, e1);
        // Convert to vector (easier to handle for now)
        std::vector<unsigned int> other(other_set.begin(), other_set.end() );

        double edge_length = (x1 - x0).norm();

        // Compute the (n-1)-dimensional covolume.
        const Eigen::Vector3d &other0 = local_node_coords[other[0]];
        const Eigen::Vector3d &other1 = local_node_coords[other[1]];
        double covolume = this->compute_covolume3d_(cc, x0, x1, other0, other1);
        // Throw an exception for 3D volumes.
        // To compute the average of the thicknesses of a control volume, one
        // has to loop over all the edges and add the thickness value to both
        // endpoints.  Then eventually, for each node, divide the resulting sum
        // by the number of connections (=number of faces of the finite
        // volume).  However, looping over edges is not (yet) possible. Hence,
        // we loop over all the cells here. This way, the edges are counted
        // several times, but it is difficult to determine how many times
        // exactly.
        //TEUCHOS_TEST_FOR_EXCEPTION(
        //    true,
        //    std::runtime_error,
        //    "Cannot calculate the average thickness in a 3D control volume yet."
        //    );

        // Compute the contributions to the finite volumes of the adjacent
        // edges.
        double pyramid_volume = 0.5 * edge_length * covolume / 3;
        // The EntityHandle (conn) is a local identifier, MOAB indices are
        // 1-based.
        cv_data[conn[e0] - 1] += pyramid_volume;
        cv_data[conn[e1] - 1] += pyramid_volume;
      }
    }
  }
}
// =============================================================================
double
mesh_tetra::
compute_covolume3d_(
    const Eigen::Vector3d &cc,
    const Eigen::Vector3d &x0,
    const Eigen::Vector3d &x1,
    const Eigen::Vector3d &other0,
    const Eigen::Vector3d &other1
    ) const
{
  double covolume = 0.0;

  // edge midpoint
  Eigen::Vector3d mp = 0.5 * (x0 + x1);

  // Compute the circumcenters of the adjacent faces.
  // This could be precomputed as well.
  Eigen::Vector3d cc_face0 = compute_triangle_circumcenter_(x0, x1, other0);
  Eigen::Vector3d cc_face1 = compute_triangle_circumcenter_(x0, x1, other1);

  // Compute the area of the quadrilateral.
  // There are some really tricky degenerate cases here, i.e., combinations
  // of when cc_face{0,1}, cc, sit outside of the tetrahedron.

  // Use the triangle (MP, local_nodes[other[0]], local_nodes[other[1]]) (in this
  // order) to gauge the orientation of the two triangles that compose the
  // quadrilateral.
  Eigen::Vector3d gauge = (other0 - mp).cross(other1 - mp);

  // Add the area of the first triangle (MP,cc_face0,cc).
  // This makes use of the right angles.
  double triangle_height0 = (mp - cc_face0).norm();
  double triangle_area0 = 0.5 * triangle_height0 * (cc_face0 - cc).norm();

  // Check if the orientation of the triangle (MP,cc_face0,cc) coincides with
  // the orientation of the gauge triangle. If yes, add the area, subtract
  // otherwise.
  Eigen::Vector3d triangle_normal0 = (cc_face0 - mp).cross(cc - mp);

  // copysign takes the absolute value of the first argument and the sign of
  // the second.
  covolume += copysign(triangle_area0, triangle_normal0.dot(gauge));

  // Add the area of the second triangle (MP,cc,cc_face1).
  // This makes use of the right angles.
  double triangle_height1 = (mp - cc_face1).norm();
  double triangle_area1 = 0.5 * triangle_height1 * (cc_face1 - cc).norm();

  // Check if the orientation of the triangle (MP,cc,cc_face1) coincides with
  // the orientation of the gauge triangle. If yes, add the area, subtract
  // otherwise.
  Eigen::Vector3d triangle_normal1 = (cc - mp).cross(cc_face1 - mp);

  // copysign takes the absolute value of the first argument and the sign of
  // the second.
  covolume += copysign(triangle_area1, triangle_normal1.dot(gauge));

  return covolume;
}
// =============================================================================
std::set<unsigned int>
mesh_tetra::
get_other_indices_(unsigned int e0, unsigned int e1) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_INEQUALITY(e0, !=, e1);
#endif
  // Get the two indices in [0,1,2,3] which are not e0, e1.
  int myint[] = {0, 1, 2, 3};
  std::set<unsigned int> a(myint, myint+4);
  a.erase(e0);
  a.erase(e1);
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(a.size(), 2);
#endif
  return a;
}
// =============================================================================
double
mesh_tetra::
get_tetrahedron_volume_(
    const Eigen::Vector3d &edge0,
    const Eigen::Vector3d &edge1,
    const Eigen::Vector3d &edge2
    ) const
{
  // Make sure the edges are not conplanar.
  double alpha = edge0.dot(edge1.cross(edge2));
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      fabs(alpha) / edge0.norm() / edge1.norm() / edge2.norm() < 1.0e-5,
      "Illegal mesh: tetrahedron too flat.\n"
      << "The following edges (with origin (0,0,0)) "
      << "seem to be conplanar:\n\n"
      << "  (0) " << edge0 << ",\n"
      << "  (1) " << edge1 << ",\n"
      << "  (2) " << edge2 << ",\n\n"
      << "Abort."
      );
  double vol = fabs(alpha) / 6.0;
  return vol;
}
// =============================================================================
Eigen::Vector3d
mesh_tetra::
compute_tetrahedron_circumcenter_(
  const std::vector<Eigen::Vector3d> &nodes
  ) const
{
  // http://www.cgafaq.info/wiki/Tetrahedron_Circumsphere
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(nodes.size(), 4);
#endif

  // Compute with respect to the first point.
  std::vector<Eigen::Vector3d> rel_nodes(3);
  for (int k = 0; k < 3; k++) {
    rel_nodes[k] = nodes[k+1] - nodes[0];
  }

  double omega = 2.0 * rel_nodes[0].dot(rel_nodes[1].cross(rel_nodes[2]));

  // don't divide by 0
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      fabs(omega) < 1.0e-10,
      "It seems that the nodes \n"
      << "\n"
      << "   " << nodes[0] << "\n"
      << "   " << nodes[1] << "\n"
      << "   " << nodes[2] << "\n"
      << "   " << nodes[3] << "\n"
      << "\n"
      << "do not form a proper tetrahedron. Abort."
      << std::endl
      );
  const double alpha = rel_nodes[0].squaredNorm() / omega;
  const double beta  = rel_nodes[1].squaredNorm() / omega;
  const double gamma = rel_nodes[2].squaredNorm() / omega;

  return nodes[0]
    + alpha * rel_nodes[1].cross(rel_nodes[2])
    + beta *  rel_nodes[2].cross(rel_nodes[0])
    + gamma * rel_nodes[0].cross(rel_nodes[1]);
}
// =============================================================================
std::vector<moab::EntityHandle>
mesh_tetra::
get_overlap_faces_() const
{
  std::vector<moab::EntityHandle> faces;

#if 0
  stk::mesh::Selector select_overlap_in_part =
    stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().universal_part())
    & (stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().locally_owned_part())
       |stk::mesh::Selector(io_broker_->bulk_data().mesh_meta_data().globally_shared_part()));

  stk::mesh::get_selected_entities(
      select_overlap_in_part,
      io_broker_->bulk_data().buckets(stk::topology::FACE_RANK),
      faces
      );
#endif
  return faces;
}
// =============================================================================
std::set<moab::EntityHandle>
mesh_tetra::
compute_boundary_nodes_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*compute_boundary_nodes_time_);
#endif
  std::set<moab::EntityHandle> _boundary_nodes;

#if 0
  // this takes a reaaaally long time
  stk::mesh::create_faces(io_broker_->bulk_data());

  auto my_faces = this->get_overlap_faces_();

  for (size_t k = 0; k < my_faces.size(); k++) {
    // if the face has one element, it's on the boundary
    if (io_broker_->bulk_data().num_elements(my_faces[k]) == 1) {
      moab::EntityHandle const * nodes =
        io_broker_->bulk_data().begin_nodes(my_faces[k]);
#ifndef NDEBUG
      TEUCHOS_ASSERT_EQUALITY(io_broker_->bulk_data().num_nodes(my_faces[k]), 3);
#endif
      _boundary_nodes.insert(nodes[0]);
      _boundary_nodes.insert(nodes[1]);
      _boundary_nodes.insert(nodes[2]);
    }
  }

#ifndef NDEBUG
  TEUCHOS_ASSERT_INEQUALITY(_boundary_nodes.size(), >, 0);
#endif
#endif
  return _boundary_nodes;
}
// =============================================================================
}  // namespace nosh
