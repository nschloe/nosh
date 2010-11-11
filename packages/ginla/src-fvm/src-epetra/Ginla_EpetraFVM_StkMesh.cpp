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
Ginla::EpetraFVM::StkMesh::
StkMesh( const Epetra_Comm                       & comm,
         const Teuchos::RCP<stk::mesh::MetaData> & metaData,
         const Teuchos::RCP<stk::mesh::BulkData> & bulkData,
         const Teuchos::RCP<VectorFieldType>     & coordinatesField
       ):
comm_( comm ),
metaData_( metaData ),
bulkData_( bulkData ),
coordinatesField_ ( coordinatesField ),
nodesMap_       ( this->createMap_( this->getOwnedNodes_()   ) ),
nodesOverlapMap_( this->createMap_( this->getOverlapNodes_() ) ),
complexMap_( this->createComplexMap_() ),
scaling_ ( Teuchos::tuple( 1.0, 1.0, 1.0 ) ),
fvmEntitiesUpToDate_( false ),
controlVolumes_( Teuchos::rcp( new Epetra_Vector( *nodesOverlapMap_ ) ) ),
coedgeEdgeRatio_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >(  3*this->getOwnedCells().size() ) ),
area_( 0.0 )
{
}
// =============================================================================
Ginla::EpetraFVM::StkMesh::
~StkMesh()
{
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

    std::vector<stk::mesh::Entity*> ownedNodes = this->getOwnedNodes_();

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
getCoedgeEdgeRatios() const
{
    if ( !fvmEntitiesUpToDate_)
        this->computeFvmEntities_();

    return coedgeEdgeRatio_;
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
Teuchos::Tuple<Point,3>
Ginla::EpetraFVM::StkMesh::
getNodeCoordinates( const stk::mesh::PairIterRelation & relation ) const
{
    Teuchos::Tuple<Point,3> localNodes;
    for ( int i=0; i<3; i++ )
    {
        double * node = stk::mesh::field_data( *coordinatesField_, *relation[i].entity() );
        localNodes[i] = Teuchos::tuple( node[0], node[1], node[2] );
    }
    return localNodes;
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
getOwnedNodes_() const
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
Ginla::EpetraFVM::StkMesh::
createMap_( const std::vector<stk::mesh::Entity*> & nodeList ) const
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
createComplexMap_() const
{
    std::vector<stk::mesh::Entity*> ownedNodes = this->getOwnedNodes_();

    // Create a map for real/imaginary out of this.
    int numDof = 2 * ownedNodes.size();
    std::vector<int> indices(numDof);
    for (int k=0; k < ownedNodes.size(); k++)
    {
        int globalNodeId = ownedNodes[k]->identifier() - 1;
        indices[2*k]   = 2*globalNodeId;
        indices[2*k+1] = 2*globalNodeId + 1;
    }
    return Teuchos::rcp(new Epetra_Map( -1, numDof, &(indices[0]), 0, comm_) );
}
// =============================================================================
// Teuchos::RCP<Epetra_Vector>
// Ginla::EpetraFVM::StkMesh::
// createStateVector( const std::vector<stk::mesh::Entity*> & ownedNodes
//                  ) const
// {
//   unsigned int neq = 2; // two components: real and imaginary part
//   Teuchos::RCP<Epetra_Vector> psi = Teuchos::rcp(new Epetra_Vector(*map));
//   for (int i=0; i < ownedNodes.size(); i++)
//   {
//       const int soln_gid = getDOF(*ownedNodes[i], 0);;
//       const int soln_lid = psi->Map().LID(soln_gid);
//       const double* psiR = stk::mesh::field_data(*stkMeshStruct->solution_field, *ownedNodes[i]);
//       (*psi)[soln_lid] = psiR[0];
//       const double* psiI = stk::mesh::field_data(*stkMeshStruct->solution_field, *ownedNodes[i]);
//       (*psi)[soln_lid+1] = psiI[0];
//   }
//   return psi;
// }
// =============================================================================
double
Ginla::EpetraFVM::StkMesh::
computeDomainArea_() const
{
  // Sum over all entries: The controlVolumes is non-overlapping.
  double norm1;
  controlVolumes_->Norm1( &norm1 );
  return norm1;
}
// =============================================================================
void
Ginla::EpetraFVM::StkMesh::
computeFvmEntities_() const
{
  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  TEUCHOS_ASSERT( !nodesOverlapMap_.is_null() );

  // Compute the volume of the (Voronoi) control cells for each point.
  if ( !controlVolumes_->Map().SameAs( *nodesOverlapMap_ ) )
      controlVolumes_->ReplaceMap( *nodesOverlapMap_ );
  controlVolumes_->PutScalar( 0.0 );

  // get owned elements
  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();

  TEUCHOS_ASSERT( !coedgeEdgeRatio_.is_null() );
  TEUCHOS_ASSERT_EQUALITY( coedgeEdgeRatio_.size(), 3*cells.size() );

  for (int k=0; k < cells.size(); k++)
  {
    stk::mesh::PairIterRelation rel = (*cells[k]).relations();

    TEST_FOR_EXCEPTION( rel.size() != 3,
                        std::runtime_error,
                        "Control volumes can only be constructed consistently with triangular elements."
                      );

    // fetch the nodal positions into 'localNodes'
    Teuchos::Tuple<Point,3> localNodes = this->getNodeCoordinates( rel );

    // compute the circumcenter of the cell
    Point cc = this->computeCircumcenter_( localNodes );

    coedgeEdgeRatio_[k] = Teuchos::ArrayRCP<double>( 3 );

    // iterate over the edges
    for ( int l=0; l<3; l++ )
    {
        // global indices of the local nodes
        int gid0 = (*rel[ l       ].entity()).identifier() - 1;
        int gid1 = (*rel[ (l+1)%3 ].entity()).identifier() - 1;

        // coordinates
        const Point & x0 = localNodes[l];
        const Point & x1 = localNodes[(l+1)%3];

        // edge midpoint
        Point mp = this->add_( 0.5, x0, 0.5, x1 );

        // sum into vector
        double val0 = this->getTriangleArea_( x0, cc, mp );
        TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->SumIntoGlobalValues( 1, &val0, &gid0 ) );

        double val1 = this->getTriangleArea_( x1, cc, mp );
        TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->SumIntoMyValues( 1, &val1, &gid1 ) );

        coedgeEdgeRatio_[k][l] = this->norm2_( this->add_( 1.0, mp, -1.0, cc ) )
                               / this->norm2_( this->add_( 1.0, x1, -1.0, x0 ) );
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
// void
// Ginla::EpetraFVM::StkMesh::
// sumInOverlapMap_( Teuchos::RCP<Epetra_Vector> x ) const
// {
//   Teuchos::RCP<Epetra_Map> nonoverlapMap =
//       Teuchos::rcp( new Epetra_Map( x->GlobalLength(), 0, comm_ ) );
//   Teuchos::RCP<Epetra_Vector> tmp =
//       Teuchos::rcp( new Epetra_Vector( *nonoverlapMap ) );
// 
//   // merge the stuff into a non-overlapping vector
//   Epetra_Export exporter( x->Map(),
//                           *nonoverlapMap
//                         );
//   tmp->Export( *x, exporter, Add );
// 
//   // map it back out to x
//   x->Import( *tmp, exporter, Insert );
// 
//   return;
// }
// =============================================================================
// Teuchos::RCP<Epetra_Map>
// Ginla::EpetraFVM::StkMesh::
// getElemsToNodesMap_() const
// {
//   // create list of nodes that need to be accessible from this process
// 
//   // Make sure that *all entries that belong to any of the elements in
//   // this core are accessible.
//   // First mark all the nodes that need to be accessible:
//   TEUCHOS_ASSERT( !elems_.is_null() );
//   TEUCHOS_ASSERT( !nodes_.is_null() );
//   int numNodes = nodes_.size();
//   Teuchos::Array<bool> mustBeAccessible( numNodes );
//   for ( int k=0; k<elems_.size(); k++ )
//       for ( int l=0; l<elems_[k].size(); l++ )
//           mustBeAccessible[ elems_[k][l] ] = true;
//   // now create the list
//   Teuchos::Array<int> entryList;
//   for ( int k=0; k<numNodes; k++ )
//       if ( mustBeAccessible[k] )
//           entryList.append( k );
// 
//   Teuchos::RCP<Epetra_Map> map =
//       Teuchos::rcp( new Epetra_Map( numNodes,
//                                     entryList.size(),
//                                     entryList.getRawPtr(),
//                                     0,
//                                     comm_
//                                   )
//                   );
// 
//   return map;
// }
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
getTriangleArea_( const Point & node0, const Point & node1, const Point & node2 ) const
{
    return 0.5 * this->norm2_( this->cross_( this->add_( 1.0, node1, -1.0, node0 ),
                                             this->add_( 1.0, node2, -1.0, node0 )
                                           )
                             );
}
// =============================================================================
Point
Ginla::EpetraFVM::StkMesh::
computeCircumcenter_( const Teuchos::Tuple<Point,3> & nodes ) const
{
  Point cc;

  double omega = 2.0 * pow( this->norm2_( this->cross_( this->add_( 1.0, nodes[0], -1.0, nodes[1] ),
                                                        this->add_( 1.0, nodes[1], -1.0, nodes[2] ) )
                                        ), 2 );

  // don't divide by 0
  TEST_FOR_EXCEPTION( fabs(omega) < 1.0e-10,
                      std::runtime_error,
                      "It seems that the nodes \n\n"
                      << "   " << nodes[0] << "\n"
                      << "   " << nodes[1] << "\n"
                      << "   " << nodes[2] << "\n"
                      << "\ndo not form a proper triangle. Abort."
                      << std::endl
                    );

  double alpha = this->dot_( this->add_( 1.0, nodes[1], -1.0, nodes[2] ), this->add_( 1.0, nodes[1], -1.0, nodes[2] ) )
               * this->dot_( this->add_( 1.0, nodes[0], -1.0, nodes[1] ), this->add_( 1.0, nodes[0], -1.0, nodes[2] ) )
               / omega;
  double beta  = this->dot_( this->add_( 1.0, nodes[2], -1.0, nodes[0] ), this->add_( 1.0, nodes[2], -1.0, nodes[0] ) )
               * this->dot_( this->add_( 1.0, nodes[1], -1.0, nodes[2] ), this->add_( 1.0, nodes[1], -1.0, nodes[0] ) )
               / omega;
  double gamma = this->dot_( this->add_( 1.0, nodes[0], -1.0, nodes[1] ), this->add_( 1.0, nodes[0], -1.0, nodes[1] ) )
               * this->dot_( this->add_( 1.0, nodes[2], -1.0, nodes[0] ), this->add_( 1.0, nodes[2], -1.0, nodes[1] ) )
               / omega;

  cc = this->add_( alpha, nodes[0], beta, nodes[1] );
  cc = this->add_( 1.0, cc, gamma, nodes[2] );

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