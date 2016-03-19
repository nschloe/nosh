#ifndef NOSH_MESHTRI_HPP
#define NOSH_MESHTRI_HPP
// =============================================================================
// includes

#include "mesh.hpp"

#include <moab/Core.hpp>

namespace nosh
{

class mesh_tri:
  public mesh
{

public:
  mesh_tri(
      const std::shared_ptr<const Teuchos::Comm<int>> & comm,
      const std::shared_ptr<moab::ParallelComm> & mcomm,
      const std::shared_ptr<moab::Core> & mb
      );

  virtual
  ~mesh_tri();

  virtual
  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  control_volumes() const
  {
    return control_volumes_;
  }

  virtual
  std::vector<edge_data>
  get_edge_data() const
  {
    return edge_data_;
  }

  virtual
  std::set<moab::EntityHandle>
  boundary_nodes() const
  {
    return boundary_nodes_;
  }

private:
  //! Compute the volume of the (Voronoi) control cells for each point.
  std::shared_ptr<Tpetra::Vector<double,int,int>>
  compute_control_volumes_() const;

  void
  compute_control_volumes_t_(
      Tpetra::Vector<double,int,int> & cv_overlap
      ) const;

  std::vector<mesh::edge_data>
  compute_edge_data_() const;

  std::set<moab::EntityHandle>
  compute_boundary_nodes_() const;

  double
  compute_covolume_(
      const Eigen::Vector3d &cc,
      const Eigen::Vector3d &x0,
      const Eigen::Vector3d &x1,
      const Eigen::Vector3d &other0
      ) const;

  unsigned int
  get_other_index_(unsigned int e0, unsigned int e1) const;

  std::set<unsigned int>
  get_other_indices_(unsigned int e0, unsigned int e1) const;

  Eigen::Vector3d
  compute_cell_circumcenter_(
      const Eigen::Vector3d &node0,
      const Eigen::Vector3d &node1,
      const Eigen::Vector3d &node2
      ) const;

  Eigen::VectorXd
  edge_coefficients_numerically_(
    const std::vector<Eigen::Vector3d> & edges
    ) const;

  double
  compute_covolume2d_(
      const Eigen::Vector3d &cc,
      const Eigen::Vector3d &x0,
      const Eigen::Vector3d &x1,
      const Eigen::Vector3d &other0
      ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> compute_edge_data_time_;
  const Teuchos::RCP<Teuchos::Time> compute_control_volumes_time_;
  const Teuchos::RCP<Teuchos::Time> compute_boundary_nodes_time_;
#endif

  std::shared_ptr<const Tpetra::Vector<double,int,int>> control_volumes_;
  const std::vector<edge_data> edge_data_;
  const std::set<moab::EntityHandle> boundary_nodes_;
};

} // namespace nosh
// =============================================================================
#endif // NOSH_MESHTRI_HPP
