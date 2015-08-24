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
// =============================================================================
// includes
#include <memory>
#include <string>
#include <vector>
#include <tuple>

#include <Teuchos_RCP.hpp>
//#ifdef NOSH_TEUCHOS_TIME_MONITOR
//#include <Teuchos_Time.hpp>
//#endif
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsGraph.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
//#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldTraits.hpp>

#include <Eigen/Dense>

// typedefs
typedef stk::mesh::Field<double, stk::mesh::Cartesian> vector_fieldType;
typedef stk::mesh::Field<double> ScalarFieldType;
//typedef stk::mesh::Field<int> IntScalarFieldType;
typedef std::tuple<stk::mesh::Entity, stk::mesh::Entity> edge;

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
      const std::shared_ptr<stk::io::StkMeshIoBroker> & broker
      );

  virtual
  ~mesh();

  double
  time() const
  {
    return time_;
  };

  void
  open_file(const std::string &output_file);

  void
  write(const double time = 0.0) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  get_vector(const std::string & field_name) const;

  std::shared_ptr<Tpetra::MultiVector<double,int,int>>
  get_multi_vector(const std::string & field_name) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  get_complex_vector(const std::string & field_name) const;

  std::vector<stk::mesh::Entity>
  get_owned_cells() const;

  std::vector<stk::mesh::Entity>
  get_overlap_edges() const;

  const std::vector<std::tuple<stk::mesh::Entity, stk::mesh::Entity>>
  my_edges() const
  {
    return edge_data_.edge_nodes;
  }

  std::vector<stk::mesh::Entity>
  owned_nodes() const
  {
    return owned_nodes_;
  }

  std::vector<stk::mesh::Entity>
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

  const vector_fieldType &
  get_node_field(const std::string & field_name) const;

  const Eigen::Vector3d
  get_node_value(
      const vector_fieldType & field,
      stk::mesh::Entity node_entity
      ) const
  {
    return Eigen::Vector3d(stk::mesh::field_data(field, node_entity));
  };

  uint64_t
  gid(const stk::mesh::Entity e) const
  {
    return io_broker_->bulk_data().identifier(e) - 1;
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

  void
  insert_complex_vector(
      const Tpetra::Vector<double,int,int> &psi,
      const std::string & field_name
      ) const;

public:
  virtual
  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  control_volumes() const = 0;

  virtual
  std::vector<double>
  edge_coefficients() const = 0;

  virtual
  std::set<stk::mesh::Entity>
  boundary_nodes() const = 0;

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
  const Teuchos::RCP<Teuchos::Time> getComplex_time_;
  const Teuchos::RCP<Teuchos::Time> multi_time_;
#endif

public:
  const std::shared_ptr<const Teuchos::Comm<int>> comm;

protected:
  // Apparently, process_output_request is not const. Make
  // the io_broker_ mutable so our write() can be const.
  const std::shared_ptr<stk::io::StkMeshIoBroker> io_broker_;

private:
  const std::vector<stk::mesh::Entity> owned_nodes_;

protected:
  const std::shared_ptr<const Tpetra::Map<int,int>> nodes_map_;
  const std::shared_ptr<const Tpetra::Map<int,int>> nodes_overlap_map_;

private:
  const std::shared_ptr<const Tpetra::Map<int,int>> complex_map_;
  const std::shared_ptr<const Tpetra::Map<int,int>> complex_overlap_map_;

protected:
  const edges_container edge_data_;

private:
  int output_channel_;

  double time_;

public:
  const std::vector<Teuchos::Tuple<int,2>> edge_gids;
  const std::vector<Teuchos::Tuple<int,4>> edge_gids_complex;

private:
  const std::vector<Teuchos::Tuple<int,2>>
  buildEdgeGids_() const;

  const std::vector<Teuchos::Tuple<int,4>>
  buildEdgeGidsComplex_() const;

  std::shared_ptr<stk::io::StkMeshIoBroker>
  read_(
      const std::string &file_name,
      const int index
      );

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  complexfield_to_vector_(
      const ScalarFieldType &real_field,
      const ScalarFieldType &imag_field
      ) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  field_to_vector_(const ScalarFieldType &field) const;

  std::shared_ptr<Tpetra::MultiVector<double,int,int>>
  field_to_vector_(
      const vector_fieldType &field,
      const int num_components
      ) const;

  std::vector<stk::mesh::Entity>
  buildOwnedNodes_(const stk::mesh::BulkData & myBulkData) const;

  std::vector<double>
  compute_edge_coefficients_() const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  buildEntitiesMap_(const std::vector<stk::mesh::Entity> &entityList) const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  buildComplexMap_(const std::vector<stk::mesh::Entity> &node_list) const;

  edges_container
  buildEdge_data_();
};
// -----------------------------------------------------------------------------

} // namespace nosh
// =============================================================================
#endif // NOSH_MESH_HPP
