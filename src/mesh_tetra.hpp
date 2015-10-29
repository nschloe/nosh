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

//#include <stk_mesh/base/Entity.hpp>
//#include <stk_mesh/base/CoordinateSystems.hpp>
//#include <stk_mesh/base/MetaData.hpp>
//#include <stk_io/StkMeshIoBroker.hpp>
//#include <stk_mesh/base/FieldTraits.hpp>

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
      const std::shared_ptr<const Teuchos::Comm<int>> & comm,
      const std::shared_ptr<stk::io::StkMeshIoBroker> & broker,
      const std::set<std::string> allocated_vector_names = {}
      );

  virtual
  ~mesh_tetra();

  virtual
  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  control_volumes() const
  {
    return control_volumes_;
  }

  virtual
  std::vector<double>
  edge_coefficients() const
  {
    return edge_coefficients_;
  }

  virtual
  std::set<stk::mesh::Entity>
  boundary_nodes() const
  {
    return boundary_nodes_;
  }

private:

  std::vector<stk::mesh::Entity>
  get_overlap_faces_() const;

  std::vector<double>
  compute_edge_coefficients_() const;

  std::set<stk::mesh::Entity>
  compute_boundary_nodes_() const;

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
    const std::vector<Eigen::Vector3d> &nodes
    ) const;

  Eigen::VectorXd
  edge_coefficients_cell_(
    const std::vector<Eigen::Vector3d> edges
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
      const std::vector<Eigen::Vector3d> &nodes
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
  const Teuchos::RCP<Teuchos::Time> compute_edge_coefficients_time_;
  const Teuchos::RCP<Teuchos::Time> compute_control_volumes_time_;
  const Teuchos::RCP<Teuchos::Time> compute_boundary_nodes_time_;
#endif

  const std::shared_ptr<const Tpetra::Vector<double,int,int>> control_volumes_;
  const std::vector<double> edge_coefficients_;
  const std::set<stk::mesh::Entity> boundary_nodes_;

};

} // namespace nosh
// =============================================================================
#endif // NOSH_MESHTETRA_HPP
