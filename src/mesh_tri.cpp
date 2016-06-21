#include "mesh_tri.hpp"

#include <set>
#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_VerboseObject.hpp>

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
  ,compute_edge_data_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tri::compute_edge_data"
        ))
  ,compute_control_volumes_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tri::compute_control_volumes"
        ))
  ,compute_boundary_vertices_time_(
      Teuchos::TimeMonitor::getNewTimer(
        "Nosh: mesh_tri::compute_boundary_vertices"
        ))
#endif
  ,control_volumes_(this->compute_control_volumes_())
  ,edge_data_(this->compute_edge_data_())
  ,boundary_surface_areas_(this->compute_boundary_surface_areas_())
{
}
// =============================================================================
mesh_tri::
~mesh_tri() = default;
// =============================================================================
std::vector<mesh::edge_data>
mesh_tri::
compute_edge_data_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*compute_edge_data_time_);
#endif

  moab::Range cells = this->mbw_->get_entities_by_dimension(0, 2);

  size_t num_cells = cells.size();

  size_t num_edges = relations_.edge_vertices.size();

  std::vector<edge_data> _edge_data(num_edges);

  // compute all coordinates
  std::vector<Eigen::Vector3d> edge_coords(num_edges);
  for (size_t k = 0; k < num_edges; k++) {
    auto tmp1 = std::get<0>(relations_.edge_vertices[k]);
    std::vector<double> coords0 = this->mbw_->get_coords({tmp1});

    tmp1 = std::get<1>(relations_.edge_vertices[k]);
    std::vector<double> coords1 = this->mbw_->get_coords({tmp1});

    edge_coords[k][0] = coords0[0] - coords1[0];
    edge_coords[k][1] = coords0[1] - coords1[1];
    edge_coords[k][2] = coords0[2] - coords1[2];

    _edge_data[k].length = edge_coords[k].norm();
  }

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

    Eigen::VectorXd edge_coeffs =
      edge_coefficients_numerically_(local_edge_coords);

    // Fill the edge coefficients into the vector.
    for (int i = 0; i < edge_coeffs.size(); i++) {
      const size_t edge_idx = this->local_index(relations_.cell_edges[k][i]);
      // const int edge_id = this->liddd
      _edge_data[edge_idx].covolume +=
        edge_coeffs[i] * _edge_data[edge_idx].length;
    }
  }

  return _edge_data;
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
  TEUCHOS_ASSERT(vertices_map_);
  TEUCHOS_ASSERT(vertices_overlap_map_);
#endif

  auto _control_volumes = std::make_shared<Tpetra::Vector<double,int,int>>(
      Teuchos::rcp(vertices_map_)
      );

  // Create temporaries to hold the overlap values for control volumes.
  Tpetra::Vector<double,int,int> cv_overlap(Teuchos::rcp(vertices_overlap_map_));

  this->compute_control_volumes_t_(cv_overlap);

  // Export control volumes to a non-overlapping map, and sum the entries.
  Teuchos::RCP<const Tpetra::Export<int,int>> exporter = Tpetra::createExport(
      Teuchos::rcp(vertices_overlap_map_),
      Teuchos::rcp(vertices_map_)
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
  for (unsigned long cell : cells) {
    const auto conn = this->mbw_->get_connectivity(cell);
    const auto splitting = this->compute_triangle_splitting_(conn);

    for (int i = 0; i < 3; i++) {
      cv_data[this->local_index(conn[i])] += splitting[i];
    }
  }

  return;
}
// =============================================================================
std::vector<double>
mesh_tri::
compute_boundary_surface_areas_() const
{
  // Store data for _all_ vertices. We actually only set the boundary ones
  // though.
  // This could be organized more efficiently with MOAB tags.
  const moab::Range vertices = this->mbw_->get_entities_by_dimension(0, 0);
  const size_t num_vertices = vertices.size();
  std::vector<double> boundary_surface_areas(num_vertices);
  // initialize to 0
  std::fill(boundary_surface_areas.begin(), boundary_surface_areas.end(), 0.0);

  for (size_t k = 0; k < this->boundary_skin_.size(); k++) {
    const size_t edge_idx = this->local_index(boundary_skin_[k]);
    const auto verts = this->mbw_->get_connectivity(boundary_skin_[k]);
    // add contributions to the verts
    boundary_surface_areas[this->local_index(verts[0])] +=
      0.5 * this->edge_data_[edge_idx].length;
    boundary_surface_areas[this->local_index(verts[1])] +=
      0.5 * this->edge_data_[edge_idx].length;
  }

  return boundary_surface_areas;
}
// =============================================================================
}  // namespace nosh
