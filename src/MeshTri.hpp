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
#ifndef NOSH_MESHTRI_HPP
#define NOSH_MESHTRI_HPP
// =============================================================================
// includes
//#include <string>
//#include <vector>
//#include <tuple>
//#include <set>
//
//#include <Teuchos_RCP.hpp>
//#include <Teuchos_DefaultComm.hpp>
//#include <Teuchos_RCPStdSharedPtrConversions.hpp>
//#ifdef NOSH_TEUCHOS_TIME_MONITOR
//#include <Teuchos_Time.hpp>
//#endif
//#include <Tpetra_Vector.hpp>
//#include <Tpetra_CrsGraph.hpp>
//
//#include <stk_mesh/base/Entity.hpp>
//#include <stk_mesh/base/CoordinateSystems.hpp>
//#include <stk_mesh/base/MetaData.hpp>
//#include <stk_io/StkMeshIoBroker.hpp>
//#include <stk_mesh/base/FieldTraits.hpp>
//
//#include <Eigen/Dense>

#include "Mesh.hpp"

//// forward declarations
//namespace stk
//{
//namespace mesh
//{
//class BulkData;
//}
//} // namespace stk
//
//// typedefs
//typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
//typedef stk::mesh::Field<double> ScalarFieldType;
//typedef stk::mesh::Field<int> IntScalarFieldType;
//typedef std::tuple<stk::mesh::Entity, stk::mesh::Entity> edge;

namespace Nosh
{

class MeshTri:
  public Mesh
{

public:
  MeshTri(
      const std::shared_ptr<const Teuchos::Comm<int>> & comm,
      const std::shared_ptr<stk::io::StkMeshIoBroker> & broker
      );

  virtual
  ~MeshTri();

  virtual
  std::shared_ptr<const Tpetra::Vector<double,int,int>>
  getControlVolumes() const
  {
    return controlVolumes_;
  }

  virtual
  std::vector<double>
  getEdgeCoefficients() const
  {
    return edgeCoefficients_;
  }

  virtual
  std::set<int>
  getBoundaryNodeGids() const
  {
    return boundaryNodeGids_;
  }

private:
  //! Compute the volume of the (Voronoi) control cells for each point.
  std::shared_ptr<Tpetra::Vector<double,int,int>>
  computeControlVolumes_() const;

  void
  computeControlVolumesT_(
      Tpetra::Vector<double,int,int> & cvOverlap
      ) const;

  std::vector<double>
  computeEdgeCoefficients_() const;

  std::set<int>
  computeBoundaryNodeGids_() const;

  double
  computeCovolume_(
      const Eigen::Vector3d &cc,
      const Eigen::Vector3d &x0,
      const Eigen::Vector3d &x1,
      const Eigen::Vector3d &other0
      ) const;

  unsigned int
  getOtherIndex_(unsigned int e0, unsigned int e1) const;

  std::set<unsigned int>
  getOtherIndices_(unsigned int e0, unsigned int e1) const;

  Eigen::Vector3d
  computeCellCircumcenter_(
      const Eigen::Vector3d &node0,
      const Eigen::Vector3d &node1,
      const Eigen::Vector3d &node2
      ) const;

  Eigen::VectorXd
  getEdgeCoefficientsNumerically_(
    const std::vector<Eigen::Vector3d> edges
    ) const;

  double
  computeCovolume2d_(
      const Eigen::Vector3d &cc,
      const Eigen::Vector3d &x0,
      const Eigen::Vector3d &x1,
      const Eigen::Vector3d &other0
      ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> computeEdgeCoefficientsTime_;
  const Teuchos::RCP<Teuchos::Time> computeControlVolumesTime_;
  const Teuchos::RCP<Teuchos::Time> computeBoundaryNodesTime_;
#endif

  const std::shared_ptr<const Tpetra::Vector<double,int,int>> controlVolumes_;
  const std::vector<double> edgeCoefficients_;
  const std::set<int> boundaryNodeGids_;

};

} // namespace Nosh
// =============================================================================
#endif // NOSH_MESHTRI_HPP
