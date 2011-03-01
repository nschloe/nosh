/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
// =============================================================================
#include "Ginla_EpetraFVM_StkMesh.h"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Export.h>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Tuple.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp> // for comm_mesh_counts
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>
#include <stk_mesh/fem/FEMInterface.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/util/UseCase_mesh.hpp>
#include <Ionit_Initializer.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
// =============================================================================
// typedefs
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
// =============================================================================
Ginla::EpetraFVM::StkMesh::
StkMesh( const Epetra_Comm                       & comm,
         const Teuchos::RCP<stk::mesh::MetaData> & metaData,
         const Teuchos::RCP<stk::mesh::BulkData> & bulkData,
         const Teuchos::RCP<VectorFieldType>     & coordinatesField
//          const Teuchos::RCP<VectorFieldType>     & thicknessField
       ):
comm_( comm ),
metaData_( metaData ),
meshData_( Teuchos::rcp( new stk::io::util::MeshData() ) ),
bulkData_( bulkData ),
coordinatesField_ ( coordinatesField ),
// thicknessField_ ( thicknessField ),
nodesMap_       ( this->createNodesMap_( this->getOwnedNodes()   ) ),
nodesOverlapMap_( this->createNodesMap_( this->getOverlapNodes() ) ),
complexMap_       ( this->createComplexMap_( this->getOwnedNodes()   ) ),
complexOverlapMap_( this->createComplexMap_( this->getOverlapNodes() ) ),
scaling_ ( Teuchos::tuple( 1.0, 1.0, 1.0 ) ),
fvmEntitiesUpToDate_( false ),
controlVolumes_( Teuchos::rcp( new Epetra_Vector( *nodesMap_ ) ) ),
coareaEdgeRatio_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( this->getOwnedCells().size() ) ),
area_( 0.0 )
{

//   meshData_->m_region->field_add( Ioss::Field( "mu",
//                                                Ioss::Field::REAL,
//                                                "scalar",
//                                                Ioss::Field::REDUCTION,
//                                                1
//                                              )
//                                 );

  return;
}
// =============================================================================
Ginla::EpetraFVM::StkMesh::
~StkMesh()
{
}
// =============================================================================
void
Ginla::EpetraFVM::StkMesh::
setOutputFile( const string & outputDir,
               const string & fileBaseName
             )
{
    // prepare the data for output
#ifdef HAVE_MPI
    const Epetra_MpiComm& mpicomm = Teuchos::dyn_cast<const Epetra_MpiComm>( comm_ );
    MPI_Comm mcomm = mpicomm.Comm();
#else
    int mcomm = 1;
#endif

    // Make sure the outputDir ends in "/".
    // Dir and filename are not concatenated properly in stk::mesh,
    std::stringstream saveOutputDir;
    saveOutputDir << outputDir << "/";

    stk::io::util::create_output_mesh( fileBaseName, // filename base
                                       "e", // extension
                                       saveOutputDir.str(), // working directoru
                                       mcomm,
                                       *bulkData_,
                                       *metaData_,
                                       *meshData_
                                     );
    return;
}
// =============================================================================
const Teuchos::RCP<stk::mesh::MetaData>
Ginla::EpetraFVM::StkMesh::
getMetaData() const
{
  return metaData_;
}
// =============================================================================
const Teuchos::RCP<stk::io::util::MeshData>
Ginla::EpetraFVM::StkMesh::
getMeshData() const
{
  return meshData_;
}
// =============================================================================
const Teuchos::RCP<stk::mesh::BulkData>
Ginla::EpetraFVM::StkMesh::
getBulkData() const
{
  return bulkData_;
}
// =============================================================================
unsigned int
Ginla::EpetraFVM::StkMesh::
getNumNodes() const
{
  return nodesMap_->NumGlobalElements();
}
// =============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::EpetraFVM::StkMesh::
getControlVolumes() const
{
  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  return controlVolumes_;
}
// =============================================================================
double
Ginla::EpetraFVM::StkMesh::
getDomainArea() const
{
  if ( !fvmEntitiesUpToDate_)
        this->computeFvmEntities_();

  return area_;
}
// =============================================================================
const Epetra_Comm &
Ginla::EpetraFVM::StkMesh::
getComm() const
{
    return comm_;
}
// =============================================================================
void
Ginla::EpetraFVM::StkMesh::
scale( const Teuchos::Tuple<double,3> & newScaling )
{
    // Prevent insanely small values.
    TEST_FOR_EXCEPTION( abs(newScaling[0]) < 1.0e-5
                        || abs(newScaling[1]) < 1.0e-5
                        || abs(newScaling[2]) < 1.0e-5,
                        std::runtime_error,
                        "Trying to scale with " << newScaling << ". This is not what you want to do."
                      );

    std::vector<stk::mesh::Entity*> ownedNodes = this->getOwnedNodes();

    // adapt the position of the nodes component by component
    for ( int i=0; i<3; i++ )
    {
        if ( newScaling[i] != scaling_[i] )
        {
           double ratio = newScaling[i] / scaling_[i];
           for ( int k=0; k<ownedNodes.size(); k++ )
           {
               double* node = stk::mesh::field_data( *coordinatesField_, *ownedNodes[k] );
               node[i] *= ratio;
           }

           // store the new scaling
           scaling_[i] = newScaling[i];

           // make sure the FVM entities get updated properly
           fvmEntitiesUpToDate_ = false;
        }
    }

    return;
}
// =============================================================================
Teuchos::Tuple<double,3>
Ginla::EpetraFVM::StkMesh::
getScaling() const
{
    return scaling_;
}
// =============================================================================
Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
Ginla::EpetraFVM::StkMesh::
getCoareaEdgeRatios() const
{
    if ( !fvmEntitiesUpToDate_)
        this->computeFvmEntities_();

    return coareaEdgeRatio_;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
Ginla::EpetraFVM::StkMesh::
getOwnedCells() const
{
  // get owned elements
  stk::mesh::Selector select_owned_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                           & stk::mesh::Selector( metaData_->locally_owned_part() );
  std::vector<stk::mesh::Entity*> cells;
  stk::mesh::get_selected_entities( select_owned_in_part,
                                    bulkData_->buckets( stk::mesh::Element ),
                                    cells
                                  );
  return cells;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
Ginla::EpetraFVM::StkMesh::
getOwnedEdges() const
{
  // get owned elements
  stk::mesh::Selector select_owned_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                           & stk::mesh::Selector( metaData_->locally_owned_part() );
  std::vector<stk::mesh::Entity*> edges;
  stk::mesh::get_selected_entities( select_owned_in_part,
                                    bulkData_->buckets( stk::mesh::Edge ),
                                    edges
                                  );
  return edges;
}
// =============================================================================
Teuchos::Array<Point>
Ginla::EpetraFVM::StkMesh::
getNodeCoordinates( const stk::mesh::PairIterRelation & relation ) const
{
    unsigned int n = relation.size();
    Teuchos::Array<Point> localNodes( n );
    for ( int i=0; i<n; i++ )
    {
        double * node = stk::mesh::field_data( *coordinatesField_, *relation[i].entity() );
        localNodes[i] = Teuchos::tuple( node[0], node[1], node[2] );
    }
    return localNodes;
}
// =============================================================================
double
Ginla::EpetraFVM::StkMesh::
getThickness( const stk::mesh::PairIterRelation & relation ) const
{
//     Teuchos::Tuple<double,3> thickness;
//     for ( int i=0; i<3; i++ )
//     {
//         double * node = stk::mesh::field_data( *coordinatesField_, *relation[i].entity() );
//         thickness[i] = node[0];
//     }
//     return thickness;
  return 0.0;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
Ginla::EpetraFVM::StkMesh::
getComplexMap() const
{
    TEUCHOS_ASSERT( !complexMap_.is_null() );
    return complexMap_;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
Ginla::EpetraFVM::StkMesh::
getOwnedNodes() const
{
    stk::mesh::Selector select_owned_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                             & stk::mesh::Selector( metaData_->locally_owned_part() );

    std::vector<stk::mesh::Entity*> ownedNodes;
    stk::mesh::get_selected_entities( select_owned_in_part,
                                      bulkData_->buckets( stk::mesh::Node ),
                                      ownedNodes
                                    );
    return ownedNodes;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
Ginla::EpetraFVM::StkMesh::
getOverlapNodes() const
{
    //  overlapnodes used for overlap map -- stored for changing coords
    stk::mesh::Selector select_overlap_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                               & (   stk::mesh::Selector( metaData_->locally_owned_part() )
                                                   | stk::mesh::Selector( metaData_->globally_shared_part() )
                                                 );

    std::vector<stk::mesh::Entity*> overlapNodes;
    stk::mesh::get_selected_entities( select_overlap_in_part,
                                      bulkData_->buckets( stk::mesh::Node ),
                                      overlapNodes
                                    );
    return overlapNodes;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
Ginla::EpetraFVM::StkMesh::
createNodesMap_( const std::vector<stk::mesh::Entity*> & nodeList ) const
{
    int numNodes = nodeList.size();
    std::vector<int> indices(numNodes);
    for (int i=0; i < numNodes; i++)
        indices[i] = nodeList[i]->identifier() - 1;

    return Teuchos::rcp(new Epetra_Map( -1, numNodes, &(indices[0]), 0, comm_) );
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
Ginla::EpetraFVM::StkMesh::
createComplexMap_( const std::vector<stk::mesh::Entity*> & nodeList ) const
{
    // Create a map for real/imaginary out of this.
    int numDof = 2 * nodeList.size();
    std::vector<int> indices(numDof);
    for (int k=0; k < nodeList.size(); k++)
    {
        int globalNodeId = nodeList[k]->identifier() - 1;
        indices[2*k]   = 2*globalNodeId;
        indices[2*k+1] = 2*globalNodeId + 1;
    }
    return Teuchos::rcp(new Epetra_Map( -1, numDof, &(indices[0]), 0, comm_) );
}
// =============================================================================
unsigned int
Ginla::EpetraFVM::StkMesh::
getCellDimension( const unsigned int numLocalNodes ) const
{
  switch ( numLocalNodes )
  {
    case 3: // triangles
        return 2;
    case 4: // tetrahedra
        return 3;
    default:
        TEST_FOR_EXCEPTION( true,
                            std::runtime_error,
                            "Control volumes can only be constructed consistently with triangular or tetrahedral elements."
                          );
  }
  return 0;
}
// =============================================================================
unsigned int
Ginla::EpetraFVM::StkMesh::
getNumEdgesPerCell( unsigned int cellDimension ) const
{
  // In n-simplices, all nodes are connected with all other nodesMap.
  // Hence, numEdges==sum_{i=1}^(numLocalNodes-1) i.
  unsigned int numEdgesPerCell = 0;
  for ( unsigned int i=1; i<cellDimension+1; i++ )
      numEdgesPerCell += i;

  return numEdgesPerCell;
}
// =============================================================================
void
Ginla::EpetraFVM::StkMesh::
computeFvmEntities_() const
{

//   std::cout << "a" << std::endl;
//   stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface( *bulkData_ );
// 
//   std::cout << "b0" << std::endl;
// //   stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface( *metaData_ );
//   std::cout << "b1" << std::endl;
// 
//   stk::mesh::PartVector empty_add_parts;
//   std::cout << "c" << std::endl;
//   stk::mesh::create_adjacent_entities( *bulkData_, empty_add_parts );
//   std::cout << "d" << std::endl;
// 
//   // count the entities for the fun of it
//   std::vector<size_t> counts ;
//   stk::mesh::comm_mesh_counts( *bulkData_ , counts );
//   std::cout << counts[0] << std::endl; // nodes
//   std::cout << counts[1] << std::endl;
//   std::cout << counts[2] << std::endl;
//   std::cout << counts[3] << std::endl; // elements
//   std::cout << counts[4] << std::endl;
//   std::cout << counts[5] << std::endl;

  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  TEUCHOS_ASSERT( !nodesOverlapMap_.is_null() );

  // Compute the volume of the (Voronoi) control cells for each point.
  if ( !controlVolumes_->Map().SameAs( *nodesMap_ ) )
      TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->ReplaceMap( *nodesMap_ ) );
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->PutScalar( 0.0 ) );

  // create a temporary to hold the overlap values
  Teuchos::RCP<Epetra_Vector> cvOverlap =
      Teuchos::rcp( new Epetra_Vector( *nodesOverlapMap_ ) );

  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  TEUCHOS_ASSERT( !coareaEdgeRatio_.is_null() );
  TEUCHOS_ASSERT_EQUALITY( coareaEdgeRatio_.size(), numCells );

  // Calculate the contributions to the finite volumes and the finite volume boundary areas cell by cell.
  for (unsigned int k=0; k < cells.size(); k++)
  {
      stk::mesh::PairIterRelation rel = (*cells[k]).relations();
      unsigned int numLocalNodes = rel.size();
      unsigned int cellDimension = this->getCellDimension( numLocalNodes );

      // Confirm that we always have the same simplices.
      TEUCHOS_ASSERT_EQUALITY( numLocalNodes, cellDimension+1 );

      // Fetch the nodal positions into 'localNodes'.
      const Teuchos::Array<Point> localNodes = this->getNodeCoordinates( rel );

      // compute the circumcenter of the cell
      Point cc;
      if ( cellDimension==2 ) // triangle
          cc = this->computeTriangleCircumcenter_( localNodes );
      else if ( cellDimension==3 ) // tetrahedron
          cc = this->computeTetrahedronCircumcenter_( localNodes );
      else
          TEST_FOR_EXCEPTION( true,
                              std::runtime_error,
                              "Control volumes can only be constructed consistently with triangular or tetrahedral elements."
                            );

      coareaEdgeRatio_[k] = Teuchos::ArrayRCP<double>( this->getNumEdgesPerCell( cellDimension ) );

      // Iterate over the edges.
      // As true edge entities are not available here, loop over all pairs of nodes.
      unsigned int edgeIndex = 0;
      for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
      {
          const Point & x0 = localNodes[e0];
          const int gid0 = (*rel[e0].entity()).identifier() - 1;
          for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
          {
              const Point & x1 = localNodes[e1];
              const int gid1 = (*rel[e1].entity()).identifier() - 1;

              // edge midpoint
              Point mp = this->add_( 0.5, x0, 0.5, x1 );

              // Get the other nodes.
              Teuchos::Tuple<unsigned int,2> other = this->getOtherIndices_( e0, e1 );

              // Compute the (n-1)-dimensional coedgeVolume.
              double coedgeVolume = 0.0;
              if ( cellDimension==2 ) // triangle
              {
                  const Point & other0 = localNodes[other[0]];
                  double coedgeLength = this->norm2_( this->add_( 1.0, mp, -1.0, cc ) );

                  // The only difficulty here is to determine whether the length of coedge
                  // is to be taken positive or negative.
                  // To this end, make sure that the order (x0, cc, mp) is of the same
                  // orientation as (x0, other0, mp).
                  Point cellNormal = this->cross_( this->add_( 1.0, other0, -1.0, x0 ),
                                                   this->add_( 1.0, mp,     -1.0, x0 )
                                                 );
                  Point ccNormal = this->cross_( this->add_( 1.0, cc, -1.0, x0 ),
                                                 this->add_( 1.0, mp, -1.0, x0 )
                                               );

                  coedgeVolume = copysign( coedgeLength, this->dot_( ccNormal, cellNormal ) );
              }
              else if ( cellDimension==3 ) // tetrahedron
              {
                  const Point & other0 = localNodes[other[0]];
                  const Point & other1 = localNodes[other[1]];

                  // Compute the circumcenters of the adjacent faces.
                  // This could be precomputed as well.
                  Point ccFace0 = this->computeTriangleCircumcenter_( x0, x1, other0 );
                  Point ccFace1 = this->computeTriangleCircumcenter_( x0, x1, other1 );

                  // Compute the area of the quadrilateral.
                  // There are some really tricky degenerate cases here, i.e., combinations
                  // of when ccFace{0,1}, cc, sit outside of the tetrahedron.

                  // Use the triangle (MP, localNodes[other[0]], localNodes[other[1]] ) (in this order)
                  // to gauge the orientation of the two triangles that compose the quadrilateral.
                  Point gauge = this->cross_( this->add_( 1.0, other0, -1.0, mp ),
                                              this->add_( 1.0, other1, -1.0, mp )
                                            );

                  // Add the area of the first triangle (MP,ccFace0,cc).
                  // This makes use of the right angles.
                  double triangleHeight0 = this->norm2_( this->add_( 1.0, mp, -1.0, ccFace0 ) );
                  double triangleArea0 = 0.5
                                      * triangleHeight0
                                      * this->norm2_( this->add_( 1.0, ccFace0, -1.0, cc ) );

                  // Check if the orientation of the triangle (MP,ccFace0,cc) coincides with
                  // the orientation of the gauge triangle. If yes, add the area, subtract otherwise.
                  Point triangleNormal0 = this->cross_( this->add_( 1.0, ccFace0, -1.0, mp ),
                                                        this->add_( 1.0, cc,      -1.0, mp )
                                                      );
                  // copysign takes the absolute value of the first argument and the sign of the second.
                  coedgeVolume += copysign( triangleArea0, this->dot_( triangleNormal0, gauge ) );

                  // Add the area of the second triangle (MP,cc,ccFace1).
                  // This makes use of the right angles.
                  double triangleHeight1 = this->norm2_( this->add_( 1.0, mp, -1.0, ccFace1 ) );
                  double triangleArea1 = 0.5
                                      * triangleHeight1
                                      * this->norm2_( this->add_( 1.0, ccFace1, -1.0, cc ) );

                  // Check if the orientation of the triangle (MP,cc,ccFace1) coincides with
                  // the orientation of the gauge triangle. If yes, add the area, subtract otherwise.
                  Point triangleNormal1 = this->cross_( this->add_( 1.0, cc,      -1.0, mp ),
                                                        this->add_( 1.0, ccFace1, -1.0, mp )
                                                      );
                  // copysign takes the absolute value of the first argument and the sign of the second.
                  coedgeVolume += copysign( triangleArea1, this->dot_( triangleNormal1, gauge ) );
              }
              else
              {
                  TEST_FOR_EXCEPTION( true,
                      std::runtime_error,
                      "Control volumes can only be constructed consistently with triangular or tetrahedral elements."
                    );
              }

              double edgeLength = this->norm2_( this->add_( 1.0, x1, -1.0, x0 ) );
              coareaEdgeRatio_[k][edgeIndex++] = coedgeVolume / edgeLength;

              // Compute the contributions to the finite volumes of the adjacent edges.
              double pyramidVolume = 0.5*edgeLength * coedgeVolume / cellDimension;
              TEUCHOS_ASSERT_EQUALITY( 0, cvOverlap->SumIntoGlobalValues( 1, &pyramidVolume, &gid0 ) );
              TEUCHOS_ASSERT_EQUALITY( 0, cvOverlap->SumIntoGlobalValues( 1, &pyramidVolume, &gid1 ) );
          }
      }
  }

  // Export control volumes to a non-overlapping map, and sum the entries.
  Epetra_Export exporter( *nodesOverlapMap_, *nodesMap_ );
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->Export( *cvOverlap, exporter, Add ) );

  // update the domain area value
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->Norm1( &area_ ) );

  fvmEntitiesUpToDate_ = true;

  return;
}
// =============================================================================
Teuchos::Tuple<unsigned int,2>
Ginla::EpetraFVM::StkMesh::
getOtherIndices_( unsigned int e0, unsigned int e1 ) const
{
  // Get the two indices in [0,1,2,3] which are not e0, e1.
  unsigned int count = 0;
  Teuchos::Tuple<unsigned int,2> otherInd;
  for ( unsigned int k=0; k<4; k++ )
  {
      if ( k!=e0 && k!=e1 )
          otherInd[count++] = k;
      TEUCHOS_ASSERT_INEQUALITY( count, <=, 2 );
  }
  return otherInd;
}
// =============================================================================
Point
Ginla::EpetraFVM::StkMesh::
add_( double alpha, const Point & x,
      double beta,  const Point & y
    ) const
{
  Point z;
  for ( int k=0; k<z.size(); k++ )
      z[k] = alpha*x[k] + beta*y[k];

  return z;
}
// =============================================================================
double
Ginla::EpetraFVM::StkMesh::
getTriangleArea_( const Point & node0,
                  const Point & node1,
                  const Point & node2
                ) const
{
    return 0.5 * this->norm2_( this->cross_( this->add_( 1.0, node1, -1.0, node0 ),
                                             this->add_( 1.0, node2, -1.0, node0 )
                                           )
                             );
}
// =============================================================================
Point
Ginla::EpetraFVM::StkMesh::
computeTriangleCircumcenter_( const Teuchos::Array<Point> & nodes ) const
{
  TEUCHOS_ASSERT_EQUALITY( nodes.size(), 3 );
  return this->computeTriangleCircumcenter_( nodes[0], nodes[1], nodes[2] );
}
// =============================================================================
Point
Ginla::EpetraFVM::StkMesh::
computeTriangleCircumcenter_( const Point & node0,
                              const Point & node1,
                              const Point & node2
                            ) const
{
  Point cc;

  double omega = 2.0 * this->norm2squared_( this->cross_( this->add_( 1.0, node0, -1.0, node1 ),
                                                          this->add_( 1.0, node1, -1.0, node2 ) )
                                          );

  // don't divide by 0
  TEST_FOR_EXCEPTION( fabs(omega) < 1.0e-10,
                      std::runtime_error,
                      "It seems that the nodes \n\n"
                      << "   " << node0 << "\n"
                      << "   " << node1 << "\n"
                      << "   " << node2 << "\n"
                      << "\ndo not form a proper triangle. Abort."
                      << std::endl
                    );

  double alpha = this->dot_( this->add_( 1.0, node1, -1.0, node2 ), this->add_( 1.0, node1, -1.0, node2 ) )
               * this->dot_( this->add_( 1.0, node0, -1.0, node1 ), this->add_( 1.0, node0, -1.0, node2 ) )
               / omega;
  double beta  = this->dot_( this->add_( 1.0, node2, -1.0, node0 ), this->add_( 1.0, node2, -1.0, node0 ) )
               * this->dot_( this->add_( 1.0, node1, -1.0, node2 ), this->add_( 1.0, node1, -1.0, node0 ) )
               / omega;
  double gamma = this->dot_( this->add_( 1.0, node0, -1.0, node1 ), this->add_( 1.0, node0, -1.0, node1 ) )
               * this->dot_( this->add_( 1.0, node2, -1.0, node0 ), this->add_( 1.0, node2, -1.0, node1 ) )
               / omega;

  cc = this->add_( alpha, node0, beta, node1 );
  cc = this->add_( 1.0, cc, gamma, node2 );

  return cc;
}
// =============================================================================
Point
Ginla::EpetraFVM::StkMesh::
computeTetrahedronCircumcenter_( const Teuchos::Array<Point> & nodes ) const
{
  // http://www.cgafaq.info/wiki/Tetrahedron_Circumsphere
  TEUCHOS_ASSERT_EQUALITY( nodes.size(), 4 );

  // Compute with respect to the first point.
  Teuchos::Array<Point> relNodes(3);
  for ( int k=0; k<3; k++ )
      relNodes[k] = this->add_( 1.0, nodes[k+1], -1.0, nodes[0] );


  double omega = 2.0 * this->dot_( relNodes[0], this->cross_( relNodes[1], relNodes[2] ) );

  // don't divide by 0
  TEST_FOR_EXCEPTION( fabs(omega) < 1.0e-10,
                      std::runtime_error,
                      "It seems that the nodes \n\n"
                      << "   " << nodes[0] << "\n"
                      << "   " << nodes[1] << "\n"
                      << "   " << nodes[2] << "\n"
                      << "   " << nodes[3] << "\n"
                      << "\ndo not form a proper tetrahedron. Abort."
                      << std::endl
                    );
  double alpha = this->norm2squared_( relNodes[0] ) / omega;
  double beta  = this->norm2squared_( relNodes[1] ) / omega;
  double gamma = this->norm2squared_( relNodes[2] ) / omega;

  Point cc;
  cc = this->add_( alpha, this->cross_( relNodes[1], relNodes[2] ),
                   beta,  this->cross_( relNodes[2], relNodes[0] ) );
  cc = this->add_( 1.0,   cc,
                   gamma, this->cross_( relNodes[0], relNodes[1] ) );

  cc = this->add_( 1.0, cc,
                   1.0, nodes[0] );

  return cc;
}
// =============================================================================
double
Ginla::EpetraFVM::StkMesh::
dot_( const Point & v, const Point & w
    ) const
{
   double sum = 0.0;
   for ( int k=0; k<v.size(); k++ )
       sum += v[k] * w[k];
   return sum;
}
// =============================================================================
Point
Ginla::EpetraFVM::StkMesh::
cross_( const Point & v, const Point & w
      ) const
{
  Point z;

  z[0] = v[1]*w[2] - v[2]*w[1];
  z[1] = v[2]*w[0] - v[0]*w[2];
  z[2] = v[0]*w[1] - v[1]*w[0];

  return z;
}
// =============================================================================
double
Ginla::EpetraFVM::StkMesh::
norm2_( const Point & x
      ) const
{
  double sum = 0.0;
  for ( int k=0; k<x.size(); k++ )
      sum += x[k]*x[k];
  return sqrt( sum );
}
// =============================================================================
double
Ginla::EpetraFVM::StkMesh::
norm2squared_( const Point & x
             ) const
{
    return this->dot_( x, x );
}
// =============================================================================
