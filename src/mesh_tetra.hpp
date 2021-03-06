#ifndef NOSH_MESHTETRA_HPP
#define NOSH_MESHTETRA_HPP
// =============================================================================
//// includes
#include <vector>
#include <set>

//#include <Teuchos_DefaultComm.hpp>
//#ifdef NOSH_TEUCHOS_TIME_MONITOR
//#include <Teuchos_Time.hpp>
//#endif
#include <Tpetra_Vector.hpp>

#include <Eigen/Dense>

#include "mesh.hpp"

#include <set>

namespace nosh
{

class mesh_tetra:
  public mesh
{

public:
  mesh_tetra(
      const std::shared_ptr<const Teuchos::Comm<int>> & _comm,
      const std::shared_ptr<moab::ParallelComm> & mcomm,
      const std::shared_ptr<moab::Core> & mb
      );

  ~mesh_tetra() override;

  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  control_volumes() const override
  {
    return control_volumes_;
  }

  std::vector<edge_data>
  get_edge_data() const override
  {
    return edge_data_;
  }

  std::vector<double>
  boundary_surface_areas() const override
  {
    return boundary_surface_areas_;
  }

private:

  std::vector<moab::EntityHandle>
  get_overlap_faces_() const;

  std::vector<mesh::edge_data>
  compute_edge_data_() const;

  std::vector<double>
  compute_boundary_surface_areas_() const;

  double
  compute_covolume_(
      const Eigen::Vector3d &cc,
      const Eigen::Vector3d &x0,
      const Eigen::Vector3d &x1,
      const Eigen::Vector3d &other0,
      const Eigen::Vector3d &other1
      ) const;

  std::set<unsigned int>
  get_other_indices_(unsigned int e0, unsigned int e1) const;

  double
  get_cell_volume_(
      const Eigen::Vector3d &edge0,
      const Eigen::Vector3d &edge1,
      const Eigen::Vector3d &edge2
      ) const;

  Eigen::Vector3d
  compute_cell_circumcenter_(
    const std::vector<Eigen::Vector3d> & vertices
    ) const;

  Eigen::VectorXd
  edge_coefficients_cell_(
    const std::vector<Eigen::Vector3d> & edges
    ) const;

  //! Compute the volume of the (Voronoi) control cells for each point.
  std::shared_ptr<Tpetra::Vector<double,int,int>>
  compute_control_volumes_() const;

  void
  compute_control_volumes_t_(Tpetra::Vector<double,int,int> & cv_overlap) const;

  double
  get_tetrahedron_volume_(
      const Eigen::Vector3d &edge0,
      const Eigen::Vector3d &edge1,
      const Eigen::Vector3d &edge2
      ) const;

  Eigen::Vector3d
  compute_tetrahedron_circumcenter_(
      const std::vector<Eigen::Vector3d> &vertices
      ) const;

  double
  compute_covolume3d_(
      const Eigen::Vector3d &cc,
      const Eigen::Vector3d &x0,
      const Eigen::Vector3d &x1,
      const Eigen::Vector3d &other0,
      const Eigen::Vector3d &other1
      ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> compute_edge_data_time_;
  const Teuchos::RCP<Teuchos::Time> compute_control_volumes_time_;
#endif

  const std::shared_ptr<const Tpetra::Vector<double,int,int>> control_volumes_;
  const std::vector<edge_data> edge_data_;
  const std::vector<double> boundary_surface_areas_;
};

} // namespace nosh
// =============================================================================
#endif // NOSH_MESHTETRA_HPP
