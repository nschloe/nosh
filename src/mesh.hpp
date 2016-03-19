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

#include "moab_wrap.hpp"

typedef std::tuple<moab::EntityHandle, moab::EntityHandle> edge;

namespace nosh
{

class mesh
{
private:
  struct entity_relations {
    //! Local edge ID -> Global node IDs.
    std::vector<edge> edge_nodes;
    //! Local cell ID -> Local edge IDs.
    std::vector<std::vector<moab::EntityHandle>> cell_edges;
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
  write(const std::string & filename) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  get_vector(const std::string & field_name) const;

  std::shared_ptr<Tpetra::Vector<double,int,int>>
  get_complex_vector(const std::string & field_name) const;

  std::shared_ptr<Tpetra::MultiVector<double,int,int>>
  get_multi_vector(const std::string & field_name) const;

  const std::vector<int>
  get_owned_gids_() const;

  const std::vector<int>
  get_overlap_gids_() const;

  const std::vector<int>
  complexify_(const std::vector<int> & ids) const;

  std::shared_ptr<Tpetra::Map<int,int>>
  get_map_(const std::vector<int> & ids) const;

  const std::vector<edge>
  my_edges() const
  {
    return relations_.edge_nodes;
  }

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

  moab::Range
  get_owned_nodes() const;

  moab::EntityHandle
  gid(const moab::EntityHandle e) const
  {
    const moab::Tag gid = this->mbw_->tag_get_handle("GLOBAL_ID");
    const auto global_ids = this->mbw_->tag_get_int_data(gid, {e});
    return global_ids[0];
  }

  // Convert a moab Entity ID to a continuous index into an array.
  size_t
  local_index(const moab::EntityHandle id) const
  {
    // MOAB EntityHandles are 1-based
    return this->mbw_->mb->id_from_handle(id) - 1;
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

  std::vector<double>
  get_data(
    const std::string & tag_name,
    const moab::Range & range
    ) const;

  Eigen::Vector3d
  get_coords(
      const moab::EntityHandle vertex
      ) const;

  const std::shared_ptr<nosh::moab_wrap> mbw_;

public:
  virtual
  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  control_volumes() const = 0;

  virtual
  std::vector<double>
  edge_coefficients() const = 0;

  virtual
  std::set<moab::EntityHandle>
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
  const Teuchos::RCP<Teuchos::Time> multi_time_;
#endif

public:
  const std::shared_ptr<const Teuchos::Comm<int>> comm;

protected:
  const std::shared_ptr<moab::ParallelComm> mcomm_;

protected:
  const std::shared_ptr<const Tpetra::Map<int,int>> nodes_map_;
  const std::shared_ptr<const Tpetra::Map<int,int>> nodes_overlap_map_;

private:
  const std::shared_ptr<const Tpetra::Map<int,int>> complex_map_;
  const std::shared_ptr<const Tpetra::Map<int,int>> complex_overlap_map_;

protected:
  const entity_relations relations_;

public:
  const std::vector<Teuchos::Tuple<int,2>> edge_lids;
  const std::vector<Teuchos::Tuple<int,4>> edge_lids_complex;
  const std::vector<Teuchos::Tuple<int,2>> edge_gids;
  const std::vector<Teuchos::Tuple<int,4>> edge_gids_complex;

private:

  const std::vector<Teuchos::Tuple<int,2>>
  build_edge_lids_() const;

  const std::vector<Teuchos::Tuple<int,4>>
  build_edge_lids_complex_() const;

  const std::vector<Teuchos::Tuple<int,2>>
  build_edge_gids_() const;

  const std::vector<Teuchos::Tuple<int,4>>
  build_edge_gids_complex_() const;

  std::shared_ptr<moab::Core>
  read_(
      const std::string &file_name,
      const int index
      );

  std::vector<double>
  compute_edge_coefficients_() const;

  std::shared_ptr<const Tpetra::Map<int,int>>
  build_map_(const std::vector<moab::EntityHandle> &entityList) const;

  entity_relations
  build_entity_relations_();
};
// -----------------------------------------------------------------------------

} // namespace nosh
// =============================================================================
#endif // NOSH_MESH_HPP
