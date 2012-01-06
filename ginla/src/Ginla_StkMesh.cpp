// @HEADER
//
//    Mesh class with compatibility to stk_mesh.
//    Copyright (C) 2010--2012  Nico Schl\"omer
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

#include <Trilinos_version.h>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Export.h>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialSpdDenseSolver.hpp>

// #include <stk_mesh/fem/FEMMetaData.hpp>
// #include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
// #include <stk_mesh/base/Comm.hpp> // for comm_mesh_counts
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>
#include <stk_mesh/base/GetEntities.hpp>

// #include <stk_io/IossBridge.hpp>
// #include <Ionit_Initializer.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif

#ifdef GINLA_TEUCHOS_TIME_MONITOR
  #include <Teuchos_TimeMonitor.hpp>
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
#ifdef GINLA_TEUCHOS_TIME_MONITOR
computeEdgeCoefficientsTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: StkMesh::computeEdgeCoefficients") ),
#endif
comm_( comm ),
metaData_( metaData ),
meshData_( Teuchos::rcp( new stk::io::MeshData() ) ),
bulkData_( bulkData ),
coordinatesField_ ( coordinatesField ),
nodesMap_       ( this->createEntitiesMap_( this->getOwnedNodes()   ) ),
nodesOverlapMap_( this->createEntitiesMap_( this->getOverlapNodes() ) ),
complexMap_       ( this->createComplexMap_( this->getOwnedNodes()   ) ),
complexOverlapMap_( this->createComplexMap_( this->getOverlapNodes() ) ),
edgesOverlapMap_( Teuchos::null ),
fvmEntitiesUpToDate_( false ),
controlVolumes_( Teuchos::rcp( new Epetra_Vector( *nodesMap_ ) ) ),
controlVolumesUpToDate_( false ),
averageThickness_( Teuchos::rcp( new Epetra_Vector( *nodesMap_ ) ) ),
edgeCoefficients_( Teuchos::ArrayRCP<double>() ),
edgeCoefficientsUpToDate_( false ),
edgeCoefficientsFallback_( Teuchos::ArrayRCP<DoubleVector>() ),
edgeCoefficientsFallbackUpToDate_( false ),
createdAdjacentEntities_( false )
{
}
// =============================================================================
StkMesh::
~StkMesh()
{
}
// =============================================================================
void
StkMesh::
openOutputChannel( const string & outputDir,
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
  if ( !controlVolumesUpToDate_ )
      this->computeControlVolumes_();
  return controlVolumes_;
}
// =============================================================================
double
StkMesh::
getDomainVolume() const
{
  if ( !controlVolumesUpToDate_)
     this->computeControlVolumes_();
  // update the domain area value
  double volume;
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->Norm1( &volume ) );

  return volume;
}
// =============================================================================
const Epetra_Comm &
StkMesh::
getComm() const
{
    return comm_;
}
// =============================================================================
Teuchos::ArrayRCP<double>
StkMesh::
getEdgeCoefficients() const
{
    if ( !edgeCoefficientsUpToDate_ )
        this->computeEdgeCoefficients_();
    return edgeCoefficients_;
}
// =============================================================================
Teuchos::ArrayRCP<DoubleVector>
StkMesh::
getEdgeCoefficientsFallback() const
{
    if ( !edgeCoefficientsFallbackUpToDate_ )
        this->computeEdgeCoefficientsFallback_();
    return edgeCoefficientsFallback_;
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
getOverlapEdges() const
{
  // Create_adjacent_entities (edges) first.
  if ( !createdAdjacentEntities_ )
  {
      stk::mesh::PartVector empty_add_parts;
      stk::mesh::create_adjacent_entities( *bulkData_, empty_add_parts );
      createdAdjacentEntities_ = true;
  }

  // get overlap edges
  stk::mesh::Selector select_overlap_in_part = stk::mesh::Selector( metaData_->universal_part() )
                                              & (   stk::mesh::Selector( metaData_->locally_owned_part() )
                                                  | stk::mesh::Selector( metaData_->globally_shared_part() )
                                                );

  std::vector<stk::mesh::Entity*> edges;
  stk::mesh::get_selected_entities( select_overlap_in_part,
                                    bulkData_->buckets( metaData_->edge_rank() ),
                                    edges
                                  );
  return edges;
}
// =============================================================================
Teuchos::ArrayRCP<DoubleVector>
StkMesh::
getNodeCoordinates( const stk::mesh::PairIterRelation & nodes ) const
{
    unsigned int n = nodes.size();
    Teuchos::ArrayRCP<DoubleVector> localNodes( n );
    for ( unsigned int i=0; i<n; i++ )
    {
        double * node = stk::mesh::field_data( *coordinatesField_, *nodes[i].entity() );
        localNodes[i] = DoubleVector(3);
        for ( unsigned int k=0; k<3; k++ )
            localNodes[i][k] = node[k];
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
    TEST_FOR_EXCEPT_MSG( true,
                         "Not yet implemented." );
    return 0.0;
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
StkMesh::
getNodesMap() const
{
#ifdef _DEBUG_
    TEUCHOS_ASSERT( !nodesMap_.is_null() );
#endif
    return nodesMap_;
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
StkMesh::
getNodesOverlapMap() const
{
#ifdef _DEBUG_
    TEUCHOS_ASSERT( !nodesOverlapMap_.is_null() );
#endif
    return nodesOverlapMap_;
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
StkMesh::
getComplexNonOverlapMap() const
{
#ifdef _DEBUG_
    TEUCHOS_ASSERT( !complexMap_.is_null() );
#endif
    return complexMap_;
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
StkMesh::
getComplexOverlapMap() const
{
#ifdef _DEBUG_
    TEUCHOS_ASSERT( !complexOverlapMap_.is_null() );
#endif
    return complexOverlapMap_;
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
StkMesh::
getEdgesOverlapMap() const
{
    if ( edgesOverlapMap_.is_null() )
      edgesOverlapMap_ = this->createEntitiesMap_( this->getOverlapEdges() );
    return edgesOverlapMap_;
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
createEntitiesMap_( const std::vector<stk::mesh::Entity*> & entityList ) const
{
    int numEntities = entityList.size();
    Teuchos::Array<int> indices(numEntities);
    for (int i=0; i < numEntities; i++)
        indices[i] = entityList[i]->identifier() - 1;

    return Teuchos::rcp(new Epetra_Map( -1, numEntities, indices.getRawPtr(), 0, comm_) );
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
        TEST_FOR_EXCEPT_MSG( true,
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
computeEdgeCoefficients_() const
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*computeEdgeCoefficientsTime_);
#endif
  
  // For edge generation and looping, refer to e-mail communication with
  // Eric Cyr and Alan Williams in April 2011.
  if ( !createdAdjacentEntities_ )
  {
      stk::mesh::PartVector empty_add_parts;
      stk::mesh::create_adjacent_entities( *bulkData_, empty_add_parts );
      createdAdjacentEntities_ = true;
  }

  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  unsigned int numEdges = this->getOverlapEdges().size();

  if ( edgeCoefficients_.size() != numEdges )
      edgeCoefficients_.resize( numEdges );

  // Calculate the edge contributions edge by edge.
  for (unsigned int k=0; k < numCells; k++)
  {
      // Get the edges.
      stk::mesh::PairIterRelation localEdges =
          cells[k]->relations( metaData_->edge_rank() );

      // Get edge coordinates.
      unsigned int numLocalEdges = localEdges.size();
      Teuchos::ArrayRCP<DoubleVector> edges( numLocalEdges );
      for ( unsigned int i=0; i<numLocalEdges; i++)
      {
          // Get the two end points.
          stk::mesh::PairIterRelation endPoints =
              (*localEdges[i].entity()).relations( metaData_->node_rank() );
#ifdef _DEBUG_
          TEUCHOS_ASSERT_EQUALITY( endPoints.size(), 2 );
#endif

          // Fetch the edge coordinates.
          const Teuchos::ArrayRCP<DoubleVector> localNodeCoords =
              this->getNodeCoordinates( endPoints );
          edges[i] = localNodeCoords[1];
          edges[i] -= localNodeCoords[0];
      }

      const DoubleVector edgeCoeffs =
          getEdgeCoefficientsNumerically_( edges );

      // Fill the edge coefficients into the vector.
      for ( unsigned int i=0; i<numLocalEdges; i++ )
      {
          const int gEdgeId = (*localEdges[i].entity()).identifier() - 1;
          const int lEdgeId = this->getEdgesOverlapMap()->LID( gEdgeId );
#ifdef _DEBUG_
          TEST_FOR_EXCEPT_MSG( lEdgeId < 0,
                               "The global index " << gEdgeId
                               << " does not seem to be present on this node." );
          TEUCHOS_ASSERT_INEQUALITY( edgeCoefficients_.size(), >, lEdgeId );
#endif
          edgeCoefficients_[lEdgeId] += edgeCoeffs[i];
      }
  }

  // stk::mesh::create_adjacent_entities doesn't work for shell elements (2D)
  // and throws a std::runtime_error. In that case, or of anything else went
  // wrong, fall back to constructing the entities via the edges of the cells.
  // This means less initialization cost, but at the costly construction of the
  // kinetic energy operator, in 2D, each each has to iterated twice (once per
  // each cell that the edge is part of). It'd be worse in the 3D.

  edgeCoefficientsUpToDate_ = true;

  return;
}
// =============================================================================lvoid
void
StkMesh::
computeEdgeCoefficientsFallback_() const
{
  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  if ( edgeCoefficientsFallback_.size() != numCells )
      edgeCoefficientsFallback_.resize( numCells );

  // Calculate the edge contributions cell by cell.
  for (unsigned int k=0; k < numCells; k++)
  {
      stk::mesh::PairIterRelation localNodes =
          cells[k]->relations( metaData_->node_rank() );
      unsigned int numLocalNodes = localNodes.size();

#ifdef _DEBUG_
      // Confirm that we always have the same simplices.
      TEUCHOS_ASSERT_EQUALITY( numLocalNodes,
                               this->getCellDimension(numLocalNodes)+1
                             );
#endif

      // Fetch the nodal positions into 'localNodes'.
      const Teuchos::ArrayRCP<const DoubleVector> localNodeCoords =
          this->getNodeCoordinates( localNodes );

      // Gather the edge coordinates.
      int numLocalEdges = numLocalNodes*(numLocalNodes-1) / 2;
      Teuchos::ArrayRCP<DoubleVector> localEdgeCoords(numLocalEdges);
      unsigned int i = 0;
      for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
      {
          for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
          {
              localEdgeCoords[i] = localNodeCoords[e1];
              localEdgeCoords[i] -= localNodeCoords[e0];
              i++;
          }
      }

      edgeCoefficientsFallback_[k] =
          this->getEdgeCoefficientsNumerically_( localEdgeCoords );
  }

  edgeCoefficientsFallbackUpToDate_ = true;

  return;
}
// =============================================================================
DoubleVector
StkMesh::
getEdgeCoefficientsNumerically_( const Teuchos::ArrayRCP<const DoubleVector> edges
                               ) const
{
    int numEdges = edges.size();

    // Build an equation system for the edge coefficients alpha_k.
    // They fulfill
    //
    //    |simplex| * <u,v> = \sum_{edges e_i} alpha_i <u,e_i> <e_i,v>
    //
    // for any pair of vectors u, v in the plane of the triangle.
    //
    double vol;
    switch ( numEdges )
    {
        case 3:
            vol = getTriangleArea_( edges[0], edges[1] );
            break;
        case 6:
            // TODO Come up with a cleaner solution here.
            try
            {
                vol = getTetrahedronVolume_( edges[0], edges[1], edges[2] );
            }
            catch(...)
            {
                // If computing the volume throws an exception, then the edges
                // chosen happened to be conplanar. Changing one of those
                // fixes this.
                vol = getTetrahedronVolume_( edges[0], edges[1], edges[3] );
            }
            break;
        default:
            TEST_FOR_EXCEPT_MSG( true,
                                 "Can only handle triangles and tetrahedra." );
    }

    Teuchos::RCP<Teuchos::SerialSymDenseMatrix<int, double> > A =
            Teuchos::rcp( new Teuchos::SerialSymDenseMatrix<int, double>(numEdges) );
    Teuchos::RCP<DoubleVector> rhs =
        Teuchos::rcp( new DoubleVector(numEdges) );
    Teuchos::RCP<DoubleVector> alpha =
        Teuchos::rcp( new DoubleVector(numEdges) );

    // Build the equation system:
    // The equation
    //
    //    |simplex| ||u||^2 = \sum_i \alpha_i <u,e_i> <e_i,u>
    //
    // has to hold for all vectors u in the plane spanned by the edges,
    // particularly by the edges themselves.
    //
    // Only fill the upper part of the Hermitian matrix.
    //
    for ( int i=0; i<numEdges; i++ )
    {
        (*rhs)(i) = vol * edges[i].dot(edges[i]);
        for ( int j=i; j<numEdges; j++ )
            (*A)(i,j) = edges[i].dot(edges[j]) * edges[j].dot(edges[i]);
    }

    // Solve the equation system for the alpha_i.
    // The system is symmetric and, if the simplex is
    // not degenerate, positive definite.
    Teuchos::SerialSpdDenseSolver<int,double> solver;
    TEUCHOS_ASSERT_EQUALITY( 0, solver.setMatrix( A ) );
    TEUCHOS_ASSERT_EQUALITY( 0, solver.setVectors( alpha, rhs ) );
    if ( solver.shouldEquilibrate() )
    {
        TEUCHOS_ASSERT_EQUALITY( 0, solver.equilibrateRHS() );
#if TRILINOS_MAJOR_MINOR_VERSION > 100803
        TEUCHOS_ASSERT_EQUALITY( 0, solver.equilibrateMatrix() );
        TEUCHOS_ASSERT_EQUALITY( 0, solver.solve() );
        TEUCHOS_ASSERT_EQUALITY( 0, solver.unequilibrateLHS() );
#else
        // A bug in Trilinos<=10.8.3 makes unequilibrateLHS() fail.
        // Work around by doing the unequilibration manually.
        // Note that this relies on the scaling being
        // s(i) = 1 / sqrt(A(i,i)).
        Teuchos::RCP<DoubleVector> diagA =
            Teuchos::rcp( new DoubleVector(numEdges) );
        for ( int k=0; k<numEdges; k++ )
            (*diagA)(k) = (*A)(k,k);
        TEUCHOS_ASSERT_EQUALITY( 0, solver.equilibrateMatrix() );
        TEUCHOS_ASSERT_EQUALITY( 0, solver.solve() );
        for ( int k=0; k<numEdges; k++ )
            (*alpha)(k) = (*alpha)(k) / sqrt( (*diagA)(k) );
#endif
    }
    else
    {
        TEUCHOS_ASSERT_EQUALITY( 0, solver.solve() );
    }

//     Teuchos::ArrayRCP<double> alphaArrayRCP( numEdges );
//     for ( int k=0; k<numEdges; k++ )
//         alphaArrayRCP[k] = (*alpha)(k,0);

    // TODO Don't explicitly copy alpha on return.
    return *alpha;
}
// =============================================================================
void
StkMesh::
computeControlVolumes_() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  TEUCHOS_ASSERT( !averageThickness_.is_null() );
  TEUCHOS_ASSERT( !nodesOverlapMap_.is_null() );
#endif

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

  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  // Calculate the contributions to the finite volumes cell by cell.
  for (unsigned int k=0; k < numCells; k++)
  {
      stk::mesh::PairIterRelation localNodes =
          cells[k]->relations( metaData_->node_rank() );
      unsigned int numLocalNodes = localNodes.size();
      unsigned int cellDimension = this->getCellDimension( numLocalNodes );

#ifdef _DEBUG_
      // Confirm that we always have the same simplices.
      TEUCHOS_ASSERT_EQUALITY( numLocalNodes, cellDimension+1 );
#endif

      // Fetch the nodal positions into 'localNodes'.
      const Teuchos::ArrayRCP<const DoubleVector> localNodeCoords =
          this->getNodeCoordinates( localNodes );

      // compute the circumcenter of the cell
      DoubleVector cc;
      if ( cellDimension==2 ) // triangle
          cc = this->computeTriangleCircumcenter_( localNodeCoords );
      else if ( cellDimension==3 ) // tetrahedron
          cc = this->computeTetrahedronCircumcenter_( localNodeCoords );
      else
          TEST_FOR_EXCEPT_MSG( true,
                               "Control volumes can only be constructed "
                               << "consistently with triangular or "
                               << "tetrahedral elements."
                             );

      // Iterate over the edges.
      // As true edge entities are not available here, loop over all pairs of local nodes.
      unsigned int edgeIndex = 0;
      for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
      {
          const DoubleVector & x0 = localNodeCoords[e0];
          const int gid0 = (*localNodes[e0].entity()).identifier() - 1;
          const int lid0 = nodesOverlapMap_->LID( gid0 );
#ifdef _DEBUG_
          TEUCHOS_ASSERT_INEQUALITY( lid0, >=, 0 );
#endif
          for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
          {
              const DoubleVector & x1 = localNodeCoords[e1];
              const int gid1 = (*localNodes[e1].entity()).identifier() - 1;
              const int lid1 = nodesOverlapMap_->LID( gid1 );
#ifdef _DEBUG_
              TEUCHOS_ASSERT_INEQUALITY( lid1, >=, 0 );
#endif

              // Get the other nodes.
              Teuchos::Tuple<unsigned int,2> other = this->getOtherIndices_( e0, e1 );

              double edgeLength = this->norm2_( this->add_( 1.0, x1, -1.0, x0 ) );

              // Compute the (n-1)-dimensional covolume.
              double covolume;
              if ( cellDimension==2 ) // triangle
              {
                  const DoubleVector & other0 = localNodeCoords[other[0]];
                  covolume = this->computeCovolume2d_( cc, x0, x1, other0 );
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
                  const DoubleVector & other0 = localNodeCoords[other[0]];
                  const DoubleVector & other1 = localNodeCoords[other[1]];
                  covolume = this->computeCovolume3d_( cc, x0, x1, other0, other1 );
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
                  TEST_FOR_EXCEPT_MSG( true,
                                       "Control volumes can only be constructed "
                                       << "consistently with triangular or "
                                       << "tetrahedral elements."
                                     );
              }

              // Compute the contributions to the finite volumes of the adjacent edges.
              double pyramidVolume = 0.5*edgeLength * covolume / cellDimension;
              (*cvOverlap)[lid0] += pyramidVolume;
              (*cvOverlap)[lid1] += pyramidVolume;
          }
      }
  }

  // Export control volumes to a non-overlapping map, and sum the entries.
  Epetra_Export exporter( *nodesOverlapMap_, *nodesMap_ );
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->Export( *cvOverlap, exporter, Add ) );

  // Export control volumes to a non-overlapping map, and sum the entries.
//   TEUCHOS_ASSERT_EQUALITY( 0, averageThickness_->Export( *atOverlap, exporter, Add ) );

  // Up until this point, averageThickness_ contains the sum of the thickness at the
  // edges. Fix this to be the average thickness.
//   for ( controlvolumes )
//       averageThickness_[k] /= number of edges;

  controlVolumesUpToDate_ = true;

  return;
}
// =============================================================================
double
StkMesh::
computeCovolume2d_( const DoubleVector & cc,
                    const DoubleVector & x0,
                    const DoubleVector & x1,
                    const DoubleVector & other0
                  ) const
{
    // edge midpoint
    DoubleVector mp = this->add_( 0.5, x0, 0.5, x1 );

    double coedgeLength = this->norm2_( this->add_( 1.0, mp, -1.0, cc ) );

    // The only difficulty here is to determine whether the length of coedge
    // is to be taken positive or negative.
    // To this end, make sure that the order (x0, cc, mp) is of the same
    // orientation as (x0, other0, mp).
    DoubleVector cellNormal = this->cross_( this->add_( 1.0, other0, -1.0, x0 ),
                                     this->add_( 1.0, mp,     -1.0, x0 )
                                   );
    DoubleVector ccNormal = this->cross_( this->add_( 1.0, cc, -1.0, x0 ),
                                   this->add_( 1.0, mp, -1.0, x0 )
                                 );

    // copysign takes the absolute value of the first argument and the sign of the second.
    return copysign( coedgeLength, this->dot_( ccNormal, cellNormal ) );
}
// =============================================================================
double
StkMesh::
computeCovolume3d_( const DoubleVector & cc,
                    const DoubleVector & x0,
                    const DoubleVector & x1,
                    const DoubleVector & other0,
                    const DoubleVector & other1
                  ) const
{
    double covolume = 0.0;

    // edge midpoint
    DoubleVector mp = this->add_( 0.5, x0, 0.5, x1 );

    // Compute the circumcenters of the adjacent faces.
    // This could be precomputed as well.
    DoubleVector ccFace0 = this->computeTriangleCircumcenter_( x0, x1, other0 );
    DoubleVector ccFace1 = this->computeTriangleCircumcenter_( x0, x1, other1 );

    // Compute the area of the quadrilateral.
    // There are some really tricky degenerate cases here, i.e., combinations
    // of when ccFace{0,1}, cc, sit outside of the tetrahedron.

    // Use the triangle (MP, localNodes[other[0]], localNodes[other[1]] ) (in this order)
    // to gauge the orientation of the two triangles that compose the quadrilateral.
    DoubleVector gauge = this->cross_( this->add_( 1.0, other0, -1.0, mp ),
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
    DoubleVector triangleNormal0 = this->cross_( this->add_( 1.0, ccFace0, -1.0, mp ),
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
    DoubleVector triangleNormal1 = this->cross_( this->add_( 1.0, cc,      -1.0, mp ),
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
#ifdef _DEBUG_
      TEUCHOS_ASSERT_INEQUALITY( count, <=, 2 );
#endif
  }
  return otherInd;
}
// =============================================================================
double
StkMesh::
getTriangleArea_( const DoubleVector & node0,
                  const DoubleVector & node1,
                  const DoubleVector & node2
                ) const
{
    return this->getTriangleArea_( this->add_( 1.0, node1, -1.0, node0 ),
                                   this->add_( 1.0, node2, -1.0, node0 )
                                 );
}
// =============================================================================
double
StkMesh::
getTriangleArea_( const DoubleVector & edge0,
                  const DoubleVector & edge1
                ) const
{
    return 0.5 * this->norm2_( this->cross_( edge0, edge1 ) );
}
// =============================================================================
double
StkMesh::
getTetrahedronVolume_( const DoubleVector & node0,
                       const DoubleVector & node1,
                       const DoubleVector & node2,
                       const DoubleVector & node3
                     ) const
{
    return this->getTetrahedronVolume_( this->add_( 1.0, node1, -1.0, node0 ),
                                        this->add_( 1.0, node2, -1.0, node0 ),
                                        this->add_( 1.0, node3, -1.0, node0 )
                                      );
}
// =============================================================================
double
StkMesh::
getTetrahedronVolume_( const DoubleVector & edge0,
                       const DoubleVector & edge1,
                       const DoubleVector & edge2
                     ) const
{
    // Make sure the edges are not conplanar.
    double alpha = edge0.dot( this->cross_( edge1, edge2 ) );
    TEST_FOR_EXCEPT_MSG( fabs(alpha) / this->norm2_(edge0)
                                     / this->norm2_(edge1)
                                     / this->norm2_(edge2)
                          < 1.0e-5,
                         "The following edges seem to be conplanar:"
                          << "\n(0) " << edge0
                          << "\n(1) " << edge1
                          << "\n(2) " << edge2 );
    double vol = fabs( alpha ) / 6.0;
    return vol;
}
// =============================================================================
DoubleVector
StkMesh::
computeTriangleCircumcenter_( const Teuchos::ArrayRCP<const DoubleVector> & nodes ) const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( nodes.size(), 3 );
#endif
  return this->computeTriangleCircumcenter_( nodes[0], nodes[1], nodes[2] );
}
// =============================================================================
DoubleVector
StkMesh::
computeTriangleCircumcenter_( const DoubleVector & node0,
                              const DoubleVector & node1,
                              const DoubleVector & node2
                            ) const
{
  DoubleVector cc;

  double omega = 2.0 * this->norm2squared_( this->cross_( this->add_( 1.0, node0, -1.0, node1 ),
                                                          this->add_( 1.0, node1, -1.0, node2 ) )
                                          );

  // don't divide by 0
  TEST_FOR_EXCEPT_MSG( fabs(omega) < 1.0e-10,
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
DoubleVector
StkMesh::
computeTetrahedronCircumcenter_( const Teuchos::ArrayRCP<const DoubleVector> & nodes ) const
{
  // http://www.cgafaq.info/wiki/Tetrahedron_Circumsphere
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( nodes.size(), 4 );
#endif

  // Compute with respect to the first point.
  Teuchos::Array<DoubleVector> relNodes(3);
  for ( int k=0; k<3; k++ )
      relNodes[k] = this->add_( 1.0, nodes[k+1], -1.0, nodes[0] );


  double omega = 2.0 * this->dot_( relNodes[0], this->cross_( relNodes[1], relNodes[2] ) );

  // don't divide by 0
  TEST_FOR_EXCEPT_MSG( fabs(omega) < 1.0e-10,
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

  DoubleVector cc;
  cc = this->add_( alpha, this->cross_( relNodes[1], relNodes[2] ),
                   beta,  this->cross_( relNodes[2], relNodes[0] ) );
  cc = this->add_( 1.0,   cc,
                   gamma, this->cross_( relNodes[0], relNodes[1] ) );

  cc = this->add_( 1.0, cc,
                   1.0, nodes[0] );

  return cc;
}
// =============================================================================
DoubleVector
StkMesh::
add_( double alpha, const DoubleVector & x,
      double beta,  const DoubleVector & y
    ) const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( x.length(), 3 );
  TEUCHOS_ASSERT_EQUALITY( y.length(), 3 );
#endif
  DoubleVector z(3);
  for ( int k=0; k<z.length(); k++ )
      z[k] = alpha*x[k] + beta*y[k];

  return z;
}
// =============================================================================
double
StkMesh::
dot_( const DoubleVector & v,
      const DoubleVector & w
    ) const
{
   double sum = 0.0;
   for ( int k=0; k<v.length(); k++ )
       sum += v[k] * w[k];
   return sum;
}
// =============================================================================
DoubleVector
StkMesh::
cross_( const DoubleVector & v,
        const DoubleVector & w
      ) const
{
  DoubleVector z(3);

  z[0] = v[1]*w[2] - v[2]*w[1];
  z[1] = v[2]*w[0] - v[0]*w[2];
  z[2] = v[0]*w[1] - v[1]*w[0];

  return z;
}
// =============================================================================
double
StkMesh::
norm2_( const DoubleVector & x
      ) const
{
  return sqrt( this->dot_(x,x) );
}
// =============================================================================
double
StkMesh::
norm2squared_( const DoubleVector & x
             ) const
{
    return this->dot_( x, x );
}
// =============================================================================
} // namespace Ginla
