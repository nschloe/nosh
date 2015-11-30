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
#ifndef NOSH_MESH_HPP
#define NOSH_MESH_HPP

// includes
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <Teuchos_RCP.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <moab/Core.hpp>
#include <moab/MOABConfig.h>
#include <moab/ParallelComm.hpp>

#include <Eigen/Dense>

typedef std::tuple<moab::EntityHandle, moab::EntityHandle> edge;

namespace nosh
{

class mesh
{
private:
  // Keep bulk_data a pointer since its copy constructor is private; this causes
  // issues when trying to copy (or initialize) Mesh_dataContainer.
  struct edges_container {
    //! Local edge ID -> Global node IDs.
    std::vector<edge> edge_nodes;
    //! Local cell ID -> Local edge IDs.
    std::vector<std::vector<int>> cell_edges;
  };

public:
  mesh(
      const std::shared_ptr<const Teuchos::Comm<int>> & comm,
      const std::shared_ptr<moab::ParallelComm> & mcomm,
      const std::shared_ptr<moab::Core> & mb
      );

  virtual
  ~mesh();

  void
  open_file(const std::string &output_file);

  void
  write(const std::string & filename) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  get_vector(const std::string & field_name) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  get_complex_vector(const std::string & field_name) const;

  std::shared_ptr<Tpetra::MultiVector<double,int,int>>
  get_multi_vector(const std::string & field_name) const;

  const std::vector<int>
  get_owned_ids_() const;

  const std::vector<int>
  get_overlap_ids_() const;

  const std::vector<int>
  complexify_(const std::vector<int> & ids) const;

  std::shared_ptr<Tpetra::Map<int,int>>
  get_map_(const std::vector<int> & ids) const;

  std::vector<moab::EntityHandle>
  get_owned_cells() const;

  std::vector<moab::EntityHandle>
  get_overlap_edges() const;

  const std::vector<edge>
  my_edges() const
  {
    return edge_data_.edge_nodes;
  }

  std::vector<moab::EntityHandle>
  owned_nodes() const
  {
    return owned_nodes_;
  }

  std::vector<moab::EntityHandle>
  get_overlap_nodes() const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  map() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(nodes_map_);
#endif
    return nodes_map_;
  }

  std::shared_ptr<const Tpetra::Map<int,int>>
  overlap_map() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(nodes_overlap_map_);
#endif
    return nodes_overlap_map_;
  }

  std::shared_ptr<const Tpetra::Map<int,int>>
  complex_map() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(complex_map_);
#endif
    return complex_map_;
  }

  std::shared_ptr<const Tpetra::Map<int,int>>
  overlap_complex_map() const
  {
#ifndef NDEBUG
    TEUCHOS_ASSERT(complex_overlap_map_);
#endif
    return complex_overlap_map_;
  }

  //const vector_fieldType &
  //get_node_field(const std::string & field_name) const;

  //const Eigen::Vector3d
  //get_node_value(
  //    const vector_fieldType & field,
  //    moab::EntityHandle node_entity
  //    ) const
  //{
  //  return Eigen::Vector3d(stk::mesh::field_data(field, node_entity));
  //};

  moab::EntityHandle
  gid(const moab::EntityHandle e) const
  {
    throw 1;
    return e;
    //return io_broker_->bulk_data().identifier(e);
  }

  moab::EntityHandle
  lid(const moab::EntityHandle e) const
  {
    throw 1;
    return e;
  }

  Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
  build_graph() const;

  Teuchos::RCP<const Tpetra::CrsGraph<int,int>>
  build_complex_graph() const;

  void
  insert_vector(
      const Tpetra::Vector<double,int,int> &x,
      const std::string & field_name
      ) const;

public:
  virtual
  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  control_volumes() const = 0;

  virtual
  std::vector<double>
  edge_coefficients() const = 0;

  //virtual
  //std::set<moab::EntityHandle>
  //boundary_nodes() const = 0;

protected:

  Eigen::Vector3d
  compute_triangle_circumcenter_(
      const std::vector<Eigen::Vector3d> &nodes
      ) const;

  Eigen::Vector3d
  compute_triangle_circumcenter_(
      const Eigen::Vector3d &node0,
      const Eigen::Vector3d &node1,
      const Eigen::Vector3d &node2
      ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> write_time_;
  const Teuchos::RCP<Teuchos::Time> multi_time_;
#endif

public:
  const std::shared_ptr<const Teuchos::Comm<int>> comm;

protected:
  const std::shared_ptr<moab::Core> mb_;
  const std::shared_ptr<moab::ParallelComm> mcomm_;

private:
  const std::vector<moab::EntityHandle> owned_nodes_;

protected:
  const std::shared_ptr<const Tpetra::Map<int,int>> nodes_map_;
  const std::shared_ptr<const Tpetra::Map<int,int>> nodes_overlap_map_;

private:
  const std::shared_ptr<const Tpetra::Map<int,int>> complex_map_;
  const std::shared_ptr<const Tpetra::Map<int,int>> complex_overlap_map_;

protected:
  const edges_container edge_data_;

public:
  const std::vector<Teuchos::Tuple<int,2>> edge_lids;
  const std::vector<Teuchos::Tuple<int,4>> edge_lids_complex;

private:

  const std::vector<Teuchos::Tuple<int,2>>
  build_edge_lids_() const;

  const std::vector<Teuchos::Tuple<int,4>>
  build_edge_lids_complex_() const;

  std::shared_ptr<moab::Core>
  read_(
      const std::string &file_name,
      const int index
      );

  //std::shared_ptr<Tpetra::Vector<double,int,int>>
  //field_to_vector_(const ScalarFieldType &field) const;

  //std::shared_ptr<Tpetra::MultiVector<double,int,int>>
  //field_to_vector_(
  //    const vector_fieldType &field,
  //    const int num_components
  //    ) const;

  std::vector<moab::EntityHandle>
  //build_owned_nodes_(const stk::mesh::BulkData & myBulkData) const;
  build_owned_nodes_() const;

  std::vector<double>
  compute_edge_coefficients_() const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  build_map_(const std::vector<moab::EntityHandle> &entityList) const;

  edges_container
  buildEdge_data_();

  std::vector<double>
  get_data_(
    const std::string & tag_name,
    const moab::Range & range
    ) const;
};
// -----------------------------------------------------------------------------

} // namespace nosh
// =============================================================================
#endif // NOSH_MESH_HPP
