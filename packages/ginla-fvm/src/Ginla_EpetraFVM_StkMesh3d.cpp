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
#include "Ginla_EpetraFVM_StkMesh3d.h"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Export.h>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Tuple.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/util/UseCase_mesh.hpp>
#include <Ionit_Initializer.h>
// =============================================================================
// typedefs
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
// =============================================================================
Ginla::EpetraFVM::StkMesh3d::
StkMesh3d( const Epetra_Comm                       & comm,
           const Teuchos::RCP<stk::mesh::MetaData> & metaData,
           const Teuchos::RCP<stk::mesh::BulkData> & bulkData,
           const Teuchos::RCP<VectorFieldType>     & coordinatesField
//            const Teuchos::RCP<VectorFieldType>     & thicknessField
         ):
comm_( comm ),
metaData_( metaData ),
meshData_( Teuchos::rcp( new stk::io::util::MeshData() ) ),
bulkData_( bulkData ),
coordinatesField_ ( coordinatesField ),
// thicknessField_ ( thicknessField ),
nodesMap_       ( this->createNodesMap_( this->getOwnedNodes()   ) ),
nodesOverlapMap_( this->createNodesMap_( this->getOverlapNodes_() ) ),
complexMap_       ( this->createComplexMap_( this->getOwnedNodes()   ) ),
complexOverlapMap_( this->createComplexMap_( this->getOverlapNodes_() ) ),
scaling_ ( Teuchos::tuple( 1.0, 1.0, 1.0 ) ),
fvmEntitiesUpToDate_( false ),
controlVolumes_( Teuchos::rcp( new Epetra_Vector( *nodesOverlapMap_ ) ) ),
coareaEdgeRatio_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >(  6*this->getOwnedCells().size() ) ),
area_( 0.0 )
{
    // prepare the data for output
#ifdef HAVE_MPI
    const Epetra_MpiComm& mpicomm = Teuchos::dyn_cast<const Epetra_MpiComm>( comm_ );
    MPI_Comm mcomm = mpicomm.Comm();
#else
    int mcomm = 1;
#endif

//   meshData_->m_region->field_add( Ioss::Field( "mu",
//                                                Ioss::Field::REAL,
//                                                "scalar",
//                                                Ioss::Field::REDUCTION,
//                                                1
//                                              )
//                                 );

    stk::io::util::create_output_mesh( "solution", // filename base
                                       "e", // extension
                                       "", // working directoru
                                       mcomm,
                                       *bulkData_,
                                       *metaData_,
                                       *meshData_
                                     );

  return;
}
// =============================================================================
Ginla::EpetraFVM::StkMesh3d::
~StkMesh3d()
{
}
// =============================================================================
const Teuchos::RCP<stk::mesh::MetaData>
Ginla::EpetraFVM::StkMesh3d::
getMetaData() const
{
  return metaData_;
}
// =============================================================================
const Teuchos::RCP<stk::io::util::MeshData>
Ginla::EpetraFVM::StkMesh3d::
getMeshData() const
{
  return meshData_;
}
// =============================================================================
const Teuchos::RCP<stk::mesh::BulkData>
Ginla::EpetraFVM::StkMesh3d::
getBulkData() const
{
  return bulkData_;
}
// =============================================================================
unsigned int
Ginla::EpetraFVM::StkMesh3d::
getNumNodes() const
{
  return nodesMap_->NumGlobalElements();
}
// =============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::EpetraFVM::StkMesh3d::
getControlVolumes() const
{
  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  return controlVolumes_;
}
// =============================================================================
double
Ginla::EpetraFVM::StkMesh3d::
getDomainArea() const
{
  return area_;
}
// =============================================================================
const Epetra_Comm &
Ginla::EpetraFVM::StkMesh3d::
getComm() const
{
    return comm_;
}
// =============================================================================
void
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
getScaling() const
{
    return scaling_;
}
// =============================================================================
Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
Ginla::EpetraFVM::StkMesh3d::
getCoareaEdgeRatios() const
{
    if ( !fvmEntitiesUpToDate_)
        this->computeFvmEntities_();

    return coareaEdgeRatio_;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
Ginla::EpetraFVM::StkMesh3d::
getOwnedCells() const
{
  // get owned elements
  stk::mesh::Selector select_owned_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                           & stk::mesh::Selector( metaData_->locally_owned_part() );
  std::vector<stk::mesh::Entity*> cells;
  stk::mesh::get_selected_entities( select_owned_in_part,
                                    bulkData_->buckets( stk::mesh::Face ),
                                    cells
                                  );
  return cells;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
getComplexMap() const
{
    TEUCHOS_ASSERT( !complexMap_.is_null() );
    return complexMap_;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
getOverlapNodes_() const
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
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
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
double
Ginla::EpetraFVM::StkMesh3d::
computeDomainArea_() const
{
  // Sum over all entries: The controlVolumes is non-overlapping.
  double norm1;
  controlVolumes_->Norm1( &norm1 );
  return norm1;
}
// =============================================================================
void
Ginla::EpetraFVM::StkMesh3d::
computeFvmEntities_() const
{
  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  TEUCHOS_ASSERT( !nodesOverlapMap_.is_null() );

  // Compute the volume of the (Voronoi) control cells for each point.
  if ( !controlVolumes_->Map().SameAs( *nodesOverlapMap_ ) )
      controlVolumes_->ReplaceMap( *nodesOverlapMap_ );
  controlVolumes_->PutScalar( 0.0 );

//   TEUCHOS_ASSERT( !coedgeEdgeRatio_.is_null() );
//   TEUCHOS_ASSERT_EQUALITY( coedgeEdgeRatio_.size(), 3*cells.size() );

  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  // Calculate the contributions to the finite volumes and the finite volume boundary areas cell by cell.
  for (int k=0; k < cells.size(); k++)
  {
      stk::mesh::PairIterRelation rel = (*cells[k]).relations();

      TEST_FOR_EXCEPTION( rel.size() != 4,
                          std::runtime_error,
                          "Control volumes can only be constructed consistently with tetrahedral elements."
                        );

      // fetch the nodal positions into 'localNodes'
      const Teuchos::Array<Point> localNodes = this->getNodeCoordinates( rel );

      // compute the circumcenter of the cell
      Point cc = this->computeTetrahedronCircumcenter_( localNodes );

      unsigned int numEdges = 6;
      coareaEdgeRatio_[k] = Teuchos::ArrayRCP<double>( numEdges );

      // Iterate over the edges.
      // As true edge entities are not available here, loop over all pairs of nodes.
      // Note that for the tetrahedron, edges connect exactly each pair of nodes.
      unsigned int numNodes = 4;
      unsigned int edgeIndex = 0;
      for ( unsigned int e0=0; e0<numNodes; e0++ )
      {
          const Point & x0 = localNodes[e0];
          int gid0 = (*rel[e0].entity()).identifier() - 1;
          for ( unsigned int e1=e0+1; e1<numNodes; e1++ )
          {
              const Point & x1 = localNodes[e1];
              int gid1 = (*rel[e1].entity()).identifier() - 1;

              double edgeLength = this->norm2_( this->add_( 1.0, x1, -1.0, x0 ) );

              // edge midpoint
              Point mp = this->add_( 0.5, x0, 0.5, x1 );

              // Get the two other nodes.
              Teuchos::Tuple<unsigned int,2> other = this->getOtherIndices_( e0, e1 );

              // Compute the circumcenters of the adjacent faces.
              // This could be precomputed as well.
              Point ccFace0 = this->computeTriangleCircumcenter_( x0, x1, localNodes[other[0]] );
              Point ccFace1 = this->computeTriangleCircumcenter_( x0, x1, localNodes[other[1]] );

              // Compute the coarea.
              double coarea = this->getQuadrilateralArea_( mp, ccFace0, cc, ccFace1 );

              coareaEdgeRatio_[k][edgeIndex++] = coarea / edgeLength;

              // Compute the contributions to the finite volumes of the adjacent edges.
              double pyramidVolume = 0.5 * edgeLength * coarea / 3.0;
              TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->SumIntoGlobalValues( 1, &pyramidVolume, &gid0 ) );
              TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->SumIntoGlobalValues( 1, &pyramidVolume, &gid1 ) );
          }
      }
  }

  // Export control volumes to a non-overlapping map, and sum the entries.
  Epetra_Export exporter( *nodesOverlapMap_, *nodesMap_ );
  controlVolumes_->Export( *controlVolumes_, exporter, Add );

  // TODO move this to another spot
  area_ = this->computeDomainArea_();

  fvmEntitiesUpToDate_ = true;

  return;
}
// =============================================================================
Teuchos::Tuple<unsigned int,2>
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
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
double
Ginla::EpetraFVM::StkMesh3d::
getQuadrilateralArea_( const Point & node0,
                       const Point & node1,
                       const Point & node2,
                       const Point & node3
                     ) const
{
    // Nodes are expected to be provided in order.

    // squared length of the diagonals
    double p = this->norm2squared_( this->add_( 1.0, node3, -1.0, node1 ) );
    double q = this->norm2squared_( this->add_( 1.0, node2, -1.0, node0 ) );

    // square length 
    double a = this->norm2squared_( node0 );
    double b = this->norm2squared_( node1 );
    double c = this->norm2squared_( node2 );
    double d = this->norm2squared_( node3 );

    // http://mathworld.wolfram.com/Quadrilateral.html
    double alpha = b+d-a-c;
    return sqrt( 4.0*p*q - alpha*alpha ) / 4.0;
}
// =============================================================================
Point
Ginla::EpetraFVM::StkMesh3d::
computeTriangleCircumcenter_( const Point & node0,
                              const Point & node1,
                              const Point & node2
                            ) const
{
  Point cc;

  double omega = 2.0 * pow( this->norm2_( this->cross_( this->add_( 1.0, node0, -1.0, node1 ),
                                                        this->add_( 1.0, node1, -1.0, node2 ) )
                                        ), 2 );

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
Ginla::EpetraFVM::StkMesh3d::
computeTriangleCircumcenter_( const Teuchos::Array<Point> & nodes ) const
{
  TEUCHOS_ASSERT_EQUALITY( nodes.size(), 3 );
  return this->computeTriangleCircumcenter_( nodes[0], nodes[1], nodes[2] );
}
// =============================================================================
Point
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
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
Ginla::EpetraFVM::StkMesh3d::
norm2squared_( const Point & x
             ) const
{
    return this->dot_( x, x );
}
// =============================================================================