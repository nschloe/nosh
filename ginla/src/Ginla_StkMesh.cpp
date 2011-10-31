// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2010, 2011  Nico Schl\"omer
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
// =============================================================================
#include "Ginla_StkMesh.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Export.h>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp> // for comm_mesh_counts
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_io/IossBridge.hpp>
#include <Ionit_Initializer.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
// =============================================================================
// typedefs
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
// =============================================================================
namespace Ginla {
// =============================================================================
StkMesh::
StkMesh( const Epetra_Comm                               & comm,
         const Teuchos::RCP<stk::mesh::fem::FEMMetaData> & metaData,
         const Teuchos::RCP<stk::mesh::BulkData>         & bulkData,
         const Teuchos::RCP<VectorFieldType>             & coordinatesField
       ):
comm_( comm ),
metaData_( metaData ),
meshData_( Teuchos::rcp( new stk::io::MeshData() ) ),
bulkData_( bulkData ),
coordinatesField_ ( coordinatesField ),
nodesMap_       ( this->createNodesMap_( this->getOwnedNodes()   ) ),
nodesOverlapMap_( this->createNodesMap_( this->getOverlapNodes() ) ),
complexMap_       ( this->createComplexMap_( this->getOwnedNodes()   ) ),
complexOverlapMap_( this->createComplexMap_( this->getOverlapNodes() ) ),
scaling_ ( Teuchos::tuple( 1.0, 1.0, 1.0 ) ),
fvmEntitiesUpToDate_( false ),
controlVolumes_( Teuchos::rcp( new Epetra_Vector( *nodesMap_ ) ) ),
averageThickness_( Teuchos::rcp( new Epetra_Vector( *nodesMap_ ) ) ),
edgeCoefficients_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >( this->getOwnedCells().size() ) ),
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
StkMesh::
~StkMesh()
{
}
// =============================================================================
void
StkMesh::
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
    std::stringstream outputFile;
    outputFile << outputDir << "/" << fileBaseName << ".e";

    stk::io::create_output_mesh( outputFile.str(),
                                 mcomm,
                                 *bulkData_,
                                 *meshData_
                               );

    stk::io::define_output_fields( *meshData_,
                                   *metaData_
                                 );

    return;
}
// =============================================================================
const Teuchos::RCP<stk::mesh::fem::FEMMetaData>
StkMesh::
getMetaData() const
{
  return metaData_;
}
// =============================================================================
const Teuchos::RCP<stk::io::MeshData>
StkMesh::
getMeshData() const
{
  return meshData_;
}
// =============================================================================
const Teuchos::RCP<stk::mesh::BulkData>
StkMesh::
getBulkData() const
{
  return bulkData_;
}
// =============================================================================
unsigned int
StkMesh::
getNumNodes() const
{
  return nodesMap_->NumGlobalElements();
}
// =============================================================================
Teuchos::RCP<const Epetra_Vector>
StkMesh::
getControlVolumes() const
{
  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  return controlVolumes_;
}
// =============================================================================
double
StkMesh::
getDomainArea() const
{
  if ( !fvmEntitiesUpToDate_)
        this->computeFvmEntities_();

  return area_;
}
// =============================================================================
const Epetra_Comm &
StkMesh::
getComm() const
{
    return comm_;
}
// =============================================================================
void
StkMesh::
scale( const Teuchos::Tuple<double,3> & newScaling )
{
    // Prevent insanely small values.
    TEUCHOS_TEST_FOR_EXCEPT_MSG( abs(newScaling[0]) < 1.0e-5
                              || abs(newScaling[1]) < 1.0e-5
                              || abs(newScaling[2]) < 1.0e-5,
                                 "Trying to scale with " << newScaling << ". This is not what you want to do."
                               );

    std::vector<stk::mesh::Entity*> ownedNodes = this->getOwnedNodes();

    // adapt the position of the nodes component by component
    for ( int i=0; i<3; i++ )
    {
        if ( newScaling[i] != scaling_[i] )
        {
           double ratio = newScaling[i] / scaling_[i];
           for ( unsigned int k=0; k<ownedNodes.size(); k++ )
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
StkMesh::
getScaling() const
{
    return scaling_;
}
// =============================================================================
Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> >
StkMesh::
getEdgeCoefficients() const
{
    if ( !fvmEntitiesUpToDate_)
        this->computeFvmEntities_();

    return edgeCoefficients_;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
StkMesh::
getOwnedCells() const
{
  // get owned elements
  stk::mesh::Selector select_owned_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                           & stk::mesh::Selector( metaData_->locally_owned_part() );
  std::vector<stk::mesh::Entity*> cells;
  stk::mesh::get_selected_entities( select_owned_in_part,
                                    bulkData_->buckets( metaData_->element_rank() ),
                                    cells
                                  );
  return cells;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
StkMesh::
getOwnedEdges() const
{
  // get owned elements
  stk::mesh::Selector select_owned_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                           & stk::mesh::Selector( metaData_->locally_owned_part() );
  std::vector<stk::mesh::Entity*> edges;
  stk::mesh::get_selected_entities( select_owned_in_part,
                                    bulkData_->buckets( metaData_->edge_rank() ),
                                    edges
                                  );
  return edges;
}
// =============================================================================
Teuchos::Array<Point>
StkMesh::
getNodeCoordinates( const stk::mesh::PairIterRelation & relation ) const
{
    unsigned int n = relation.size();
    Teuchos::Array<Point> localNodes( n );
    for ( unsigned int i=0; i<n; i++ )
    {
        double * node = stk::mesh::field_data( *coordinatesField_, *relation[i].entity() );
        localNodes[i] = Teuchos::tuple( node[0], node[1], node[2] );
    }
    return localNodes;
}
// =============================================================================
double
StkMesh::
getThickness( const stk::mesh::PairIterRelation & relation ) const
{
//     Teuchos::Tuple<double,3> thickness;
//     for ( int i=0; i<3; i++ )
//     {
//         double * node = stk::mesh::field_data( *coordinatesField_, *relation[i].entity() );
//         thickness[i] = node[0];
//     }
//     return thickness;
    TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                                 "Not yet implemented." );
    return 0.0;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
StkMesh::
getNodesMap() const
{
    TEUCHOS_ASSERT( !nodesMap_.is_null() );
    return nodesMap_;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
StkMesh::
getNodesOverlapMap() const
{
    TEUCHOS_ASSERT( !nodesOverlapMap_.is_null() );
    return nodesOverlapMap_;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
StkMesh::
getComplexMap() const
{
    TEUCHOS_ASSERT( !complexMap_.is_null() );
    return complexMap_;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
StkMesh::
getComplexOverlapMap() const
{
    TEUCHOS_ASSERT( !complexOverlapMap_.is_null() );
    return complexOverlapMap_;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
StkMesh::
getOwnedNodes() const
{
    stk::mesh::Selector select_owned_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                             & stk::mesh::Selector( metaData_->locally_owned_part() );

    std::vector<stk::mesh::Entity*> ownedNodes;
    stk::mesh::get_selected_entities( select_owned_in_part,
                                      bulkData_->buckets( metaData_->node_rank() ),
                                      ownedNodes
                                    );

    return ownedNodes;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
StkMesh::
getOverlapNodes() const
{
    //  overlapnodes used for overlap map -- stored for changing coords
    stk::mesh::Selector select_overlap_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                               & (   stk::mesh::Selector( metaData_->locally_owned_part() )
                                                   | stk::mesh::Selector( metaData_->globally_shared_part() )
                                                 );

    std::vector<stk::mesh::Entity*> overlapNodes;
    stk::mesh::get_selected_entities( select_overlap_in_part,
                                      bulkData_->buckets( metaData_->node_rank() ),
                                      overlapNodes
                                    );

    return overlapNodes;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
StkMesh::
createNodesMap_( const std::vector<stk::mesh::Entity*> & nodeList ) const
{
    int numNodes = nodeList.size();
    Teuchos::Array<int> indices(numNodes);
    for (int i=0; i < numNodes; i++)
        indices[i] = nodeList[i]->identifier() - 1;

    return Teuchos::rcp(new Epetra_Map( -1, numNodes, indices.getRawPtr(), 0, comm_) );
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
StkMesh::
createComplexMap_( const std::vector<stk::mesh::Entity*> & nodeList ) const
{
    // Create a map for real/imaginary out of this.
    int numDof = 2 * nodeList.size();
    Teuchos::Array<int> indices(numDof);
    for ( unsigned int k=0; k < nodeList.size(); k++ )
    {
        int globalNodeId = nodeList[k]->identifier() - 1;
        indices[2*k]   = 2*globalNodeId;
        indices[2*k+1] = 2*globalNodeId + 1;
    }
    return Teuchos::rcp(new Epetra_Map( -1, numDof, indices.getRawPtr(), 0, comm_) );
}
// =============================================================================
unsigned int
StkMesh::
getCellDimension( const unsigned int numLocalNodes ) const
{
  switch ( numLocalNodes )
  {
    case 3: // triangles
        return 2;
    case 4: // tetrahedra
        return 3;
    default:
        TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                                     "Control volumes can only be constructed "
                                     << "consistently with triangular or "
                                     << "tetrahedral elements."
                                   );
  }
  return 0;
}
// =============================================================================
unsigned int
StkMesh::
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
StkMesh::
computeFvmEntities_() const
{
//
// For edge generation an looping, refer to e-mail communication with Eric Cyr and Alan Williams
// in April 2011.
//
//   stk::mesh::PartVector empty_add_parts;
//   stk::mesh::create_adjacent_entities( *bulkData_, empty_add_parts );
//
//   std::vector<stk::mesh::Entity*> edges = this->getOwnedEdges();
//
//   std::cout << "My edges: " << edges.size() << std::endl;


  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  TEUCHOS_ASSERT( !averageThickness_.is_null() );
  TEUCHOS_ASSERT( !nodesOverlapMap_.is_null() );

  // Compute the volume of the (Voronoi) control cells for each point.
  if ( !controlVolumes_->Map().SameAs( *nodesMap_ ) )
      TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->ReplaceMap( *nodesMap_ ) );
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->PutScalar( 0.0 ) );

  // Compute the average thickness for each control volume.
  if ( !averageThickness_->Map().SameAs( *nodesMap_ ) )
      TEUCHOS_ASSERT_EQUALITY( 0, averageThickness_->ReplaceMap( *nodesMap_ ) );
  TEUCHOS_ASSERT_EQUALITY( 0, averageThickness_->PutScalar( 0.0 ) );

  // Create temporaries to hold the overlap values for control volumes and
  // average thickness.
  Teuchos::RCP<Epetra_Vector> cvOverlap =
      Teuchos::rcp( new Epetra_Vector( *nodesOverlapMap_ ) );
//   Teuchos::RCP<Epetra_Vector> atOverlap =
//       Teuchos::rcp( new Epetra_Vector( *nodesOverlapMap_ ) );

  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  TEUCHOS_ASSERT( !edgeCoefficients_.is_null() );
  TEUCHOS_ASSERT_EQUALITY( edgeCoefficients_.size(), numCells );

  // Calculate the contributions to the finite volumes and the finite volume boundary areas cell by cell.
  for (unsigned int k=0; k < numCells; k++)
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
          TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                                       "Control volumes can only be constructed "
                                       << "consistently with triangular or "
                                       << "tetrahedral elements."
                                     );

//      edgeCoefficients_[k] = Teuchos::ArrayRCP<double>( this->getNumEdgesPerCell( cellDimension ) );

      edgeCoefficients_[k] = getEdgeCoefficientsNumerically_( localNodes );

      // Iterate over the edges.
      // As true edge entities are not available here, loop over all pairs of local nodes.
      unsigned int edgeIndex = 0;
      for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
      {
          const Point & x0 = localNodes[e0];
          const int gid0 = (*rel[e0].entity()).identifier() - 1;
//          const int lid0 = nodesMap_->LID( gid0 );
          for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
          {
              const Point & x1 = localNodes[e1];
              const int gid1 = (*rel[e1].entity()).identifier() - 1;
//              const int lid1 = nodesMap_->LID( gid1 );

              // Get the other nodes.
              Teuchos::Tuple<unsigned int,2> other = this->getOtherIndices_( e0, e1 );

              double edgeLength = this->norm2_( this->add_( 1.0, x1, -1.0, x0 ) );

              // Compute the (n-1)-dimensional covolume.
              double covolume;
              if ( cellDimension==2 ) // triangle
              {
                  const Point & other0 = localNodes[other[0]];
                  covolume = this->computeCovolume2d_( cc, x0, x1, other0 );

//                  std::cout << edgeCoefficients_[k][edgeIndex]
//                            << " " << covolume/edgeLength
//                            << " "  << fabs(edgeCoefficients_[k][edgeIndex]-covolume/edgeLength)
//                            << std::endl;
                  // For testing:
                  // Check if the computed coefficients actually coincide with
                  // what is computed analytically.
                  TEUCHOS_ASSERT_INEQUALITY( fabs(edgeCoefficients_[k][edgeIndex]-covolume/edgeLength),
                                             <,
                                             1.0e-12 );
                  edgeCoefficients_[k][edgeIndex++] = covolume / edgeLength;

                  // The problem with counting the average thickness in 2D is the following.
                  // Ideally, one would want to loop over all edges, add the midpoint value
                  // of the thickness to both of the edge end points, and eventually loop over
                  // all endpoints and divide by the number of edges (connections) they have
                  // with neighboring nodes).
                  // Unfortunately, this is impossible now b/c there's no edge generation
                  // for shells in Trilinos yet (2011-04-15).
                  // As a workaround, one could loop over all cells, and then all pairs of
                  // nodes to retrieve the edges. In 2D, almost all of the edges would be
                  // counted twice this way as they belong to two cells. This is true for
                  // all but the boundary edges. Again, it is difficult (impossible?) to
                  // know what the boundary edges are, and hence which values to divide by
                  // 2. Dividing them all by two would result in an artificially lower
                  // thickness near the boundaries. This is not what we want.
              }
              else if ( cellDimension==3 ) // tetrahedron
              {
                  const Point & other0 = localNodes[other[0]];
                  const Point & other1 = localNodes[other[1]];
                  covolume = this->computeCovolume3d_( cc, x0, x1, other0, other1 );


//                  // The 3d-calculation of the edge coefficients as covolume-
//                  // edge length ratio is wrong.
//                  std::cout << edgeCoefficients_[k][edgeIndex] << " " << covolume/edgeLength << std::endl;
//                  TEUCHOS_ASSERT_INEQUALITY( fabs(edgeCoefficients_[k][edgeIndex]-covolume/edgeLength),
//                                             <,
//                                             1.0e-12 );
//                  edgeCoefficients_[k][edgeIndex++] = covolume / edgeLength;

                  // Throw an exception for 3D volumes.
                  // To compute the average of the thicknesses of a control volume, one has to loop
                  // over all the edges and add the thickness value to both endpoints.
                  // Then eventually, for each node, divide the resulting sum by the number of connections
                  // (=number of faces of the finite volume).
                  // However, looping over edges is not (yet) possible. Hence, we loop over all the
                  // cells here. This way, the edges are counted several times, but it is difficult
                  // to determine how many times exactly.
//                   TEST_FOR_EXCEPTION( true,
//                                       std::runtime_error,
//                                       "Cannot calculate the average thickness in a 3D control volume yet."
//                                     );
              }
              else
              {
                  TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                                               "Control volumes can only be constructed "
                                               << "consistently with triangular or "
                                               << "tetrahedral elements."
                                             );
              }

              // Compute the contributions to the finite volumes of the adjacent edges.
              double pyramidVolume = 0.5*edgeLength * covolume / cellDimension;
              TEUCHOS_ASSERT_EQUALITY( 0, cvOverlap->SumIntoGlobalValues( 1, &pyramidVolume, &gid0 ) );
              TEUCHOS_ASSERT_EQUALITY( 0, cvOverlap->SumIntoGlobalValues( 1, &pyramidVolume, &gid1 ) );
          }
      }
  }

  // Export control volumes to a non-overlapping map, and sum the entries.
  Epetra_Export exporter( *nodesOverlapMap_, *nodesMap_ );
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->Export( *cvOverlap, exporter, Add ) );
//   TEUCHOS_ASSERT_EQUALITY( 0, averageThickness_->Export( *atOverlap, exporter, Add ) );

  // Up until this point, averageThickness_ contains the sum of the thickness at the
  // edges. Fix this to be the average thickness.
//   for ( controlvolumes )
//       averageThickness_[k] /= number of edges;

  // update the domain area value
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->Norm1( &area_ ) );

  fvmEntitiesUpToDate_ = true;

  return;
}
// =============================================================================
Teuchos::ArrayRCP<double>
StkMesh::
getEdgeCoefficientsNumerically_( const Teuchos::Array<Point> localNodes ) const
{
    // Build an equation system for the edge coefficients alpha_k.
    // They fulfill
    //
    //    |simplex| * <u,v> = \sum_{edges e_i} alpha_i <u,e_i> <e_i,v>
    //
    // for any pair of vectors u, v in the plane of the triangle.
    //
    double vol;
    int numEdges;
    unsigned int numLocalNodes = localNodes.size();
    if (localNodes.size()==3)
    {
        vol = getTriangleArea_( localNodes[0], localNodes[1], localNodes[2] );
        numEdges = 3;
    }
    else if (localNodes.size()==4)
    {
        vol = getTetrahedronVolume_( localNodes[0], localNodes[1], localNodes[2], localNodes[3] );
        numEdges = 6;
    }
    else
        TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                                     "Can only handle triangles and tetrahedra." );

//    std::cout << "\nTriangle is \n";
//    for ( int i=0; i<localNodes.size(); i++ )
//    {
//        for ( int j=0; j<3; j++ )
//            std::cout << " " << localNodes[i][j];
//        std::cout << std::endl;
//    }

    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > A =
            Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, double>(numEdges,numEdges) );
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > rhs =
        Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, double>(numEdges,1) );
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > alpha =
        Teuchos::rcp( new Teuchos::SerialDenseMatrix<int, double>(numEdges,1) );

    // Gather the edges.
    Teuchos::ArrayRCP<Teuchos::SerialDenseVector<int, double> > edges(numEdges);
    unsigned int i = 0;
    for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
    {
        for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
        {
            edges[i] = Teuchos::SerialDenseVector<int, double>(3);
            for ( int k=0; k<3; k++ )
                edges[i][k] = localNodes[e1][k] - localNodes[e0][k];
            i++;
        }
    }

    // Build the equation system:
    // The equation
    //
    //    |simplex| ||u||^2 = \sum_i \alpha_i <u,e_i> <e_i,u>
    //
    // has to hold for all vectors u in the plane spanned by the edges,
    // particularly by the edges themselves.
    //
    for ( int i=0; i<numEdges; i++ )
    {
        (*rhs)(i,0) = vol * edges[i].dot(edges[i]);
        for ( int j=0; j<numEdges; j++ )
            (*A)(i,j) = edges[i].dot(edges[j]) * edges[j].dot(edges[i]);
    }

//    A->print( std::cout );
//    rhs->print( std::cout );

    // solve the equation system for the alpha_i
    Teuchos::SerialDenseSolver<int,double> solver;
    TEUCHOS_ASSERT_EQUALITY( 0, solver.setMatrix( A ) );
    TEUCHOS_ASSERT_EQUALITY( 0, solver.setVectors( alpha, rhs ) );
    if ( false ) // ( solver.shouldEquilibrate() ) // reactivate when equilibration works in trilinos
    {
//        std::cout << "no equi" << std::endl;
//        TEUCHOS_ASSERT_EQUALITY( 0, solver.computeEquilibrateScaling() );
        TEUCHOS_ASSERT_EQUALITY( 0, solver.equilibrateMatrix() );
        TEUCHOS_ASSERT_EQUALITY( 0, solver.equilibrateRHS() );
        TEUCHOS_ASSERT_EQUALITY( 0, solver.solve() );
        TEUCHOS_ASSERT_EQUALITY( 0, solver.unequilibrateLHS() );
    }
    else
    {
//        std::cout << "no equi" << std::endl;
        TEUCHOS_ASSERT_EQUALITY( 0, solver.solve() );
    }


    Teuchos::ArrayRCP<double> alphaArrayRCP( numEdges );
    for ( int k=0; k<numEdges; k++ )
        alphaArrayRCP[k] = (*alpha)(k,0);

//    std::cout << "Got coefficients";
//    for ( int k=0; k<numEdges; k++ )
//        std::cout << " " << alphaArrayRCP[k];
//    std::cout << std::endl;

    return alphaArrayRCP;
}
// =============================================================================
double
StkMesh::
computeCovolume2d_( const Point & cc,
                    const Point & x0,
                    const Point & x1,
                    const Point & other0
                  ) const
{
    // edge midpoint
    Point mp = this->add_( 0.5, x0, 0.5, x1 );

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

    // copysign takes the absolute value of the first argument and the sign of the second.
    return copysign( coedgeLength, this->dot_( ccNormal, cellNormal ) );
}
// =============================================================================
double
StkMesh::
computeCovolume3d_( const Point & cc,
                    const Point & x0,
                    const Point & x1,
                    const Point & other0,
                    const Point & other1
                  ) const
{
    double covolume = 0.0;

    // edge midpoint
    Point mp = this->add_( 0.5, x0, 0.5, x1 );

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
    covolume += copysign( triangleArea0, this->dot_( triangleNormal0, gauge ) );

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
    covolume += copysign( triangleArea1, this->dot_( triangleNormal1, gauge ) );

    return covolume;
}
// =============================================================================
Teuchos::Tuple<unsigned int,2>
StkMesh::
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
StkMesh::
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
StkMesh::
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
StkMesh::
getTetrahedronVolume_( const Point & node0,
                       const Point & node1,
                       const Point & node2,
                       const Point & node3
                     ) const
{
    return fabs( this->dot_( this->add_( 1.0, node1, -1.0, node0 ),
                             this->cross_( this->add_( 1.0, node2, -1.0, node0 ),
                                           this->add_( 1.0, node3, -1.0, node0 ) )
                           )
               ) / 6.0;
}
// =============================================================================
Point
StkMesh::
computeTriangleCircumcenter_( const Teuchos::Array<Point> & nodes ) const
{
  TEUCHOS_ASSERT_EQUALITY( nodes.size(), 3 );
  return this->computeTriangleCircumcenter_( nodes[0], nodes[1], nodes[2] );
}
// =============================================================================
Point
StkMesh::
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
  TEUCHOS_TEST_FOR_EXCEPT_MSG( fabs(omega) < 1.0e-10,
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
StkMesh::
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
  TEUCHOS_TEST_FOR_EXCEPT_MSG( fabs(omega) < 1.0e-10,
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
StkMesh::
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
StkMesh::
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
StkMesh::
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
StkMesh::
norm2squared_( const Point & x
             ) const
{
    return this->dot_( x, x );
}
// =============================================================================
} // namespace Ginla
