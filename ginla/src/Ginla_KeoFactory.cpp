// @HEADER
//
//    Factory class that hosts the kinetic energy operator.
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
// includes
#include "Ginla_KeoFactory.hpp"
#include "Ginla_StkMesh.hpp"
#include "Ginla_Helpers.hpp"
#include "Ginla_MagneticVectorPotential_Virtual.hpp"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>

#include <Epetra_FECrsGraph.h>

#include <ml_epetra_preconditioner.h>

#ifdef GINLA_TEUCHOS_TIME_MONITOR
  #include <Teuchos_TimeMonitor.hpp>
#endif

namespace Ginla {
// =============================================================================
KeoFactory::
KeoFactory( const Teuchos::RCP<const Ginla::StkMesh>           & mesh,
            const Teuchos::RCP<const Epetra_Vector>            & thickness,
            const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp
          ):
#ifdef GINLA_TEUCHOS_TIME_MONITOR
        keoFillTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::fillKeo_") ),
#endif
        mesh_ ( mesh ),
        thickness_( thickness ),
        mvp_( mvp ),
        globalIndexCache_( Teuchos::ArrayRCP<Epetra_IntSerialDenseVector>() ),
        globalIndexCacheUpToDate_( false ),
        globalIndexFallbackCache_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<Epetra_IntSerialDenseVector> >() ),
        globalIndexFallbackCacheUpToDate_( false ),
        keoGraph_( this->buildKeoGraph_() ), // build the graph immediately
        keo_( Teuchos::rcp( new Epetra_FECrsMatrix( Copy, *keoGraph_ ) ) ),
        keoBuildParameters_( Teuchos::null ),
        keoDMu_( Teuchos::rcp( new Epetra_FECrsMatrix( Copy, *keoGraph_ ) ) ),
        keoDMuBuildParameters_( Teuchos::null ),
        keoDTheta_( Teuchos::rcp( new Epetra_FECrsMatrix( Copy, *keoGraph_ ) ) ),
        keoDThetaBuildParameters_( Teuchos::null ),
        alphaCache_( Teuchos::ArrayRCP<double>() ),
        alphaCacheUpToDate_( false ),
        alphaFallbackCache_( Teuchos::ArrayRCP<DoubleVector>() ),
        alphaFallbackCacheUpToDate_( false )
{
}
// =============================================================================
KeoFactory::
~KeoFactory()
{
}
// =============================================================================
const Epetra_Comm &
KeoFactory::
getComm() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
    return mesh_->getComm();
}
// =============================================================================
void
KeoFactory::
updateParameters( const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams
                ) const
{
  // set the parameters
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mvpParams.is_null() );
  TEUCHOS_ASSERT( !mvp_.is_null() );
#endif
  mvp_->setParameters( *mvpParams );
  return;
}
// =============================================================================
const Teuchos::RCP<const LOCA::ParameterVector>
KeoFactory::
getMvpParameters() const
{
  return mvp_->getParameters();
}
// =============================================================================
Teuchos::RCP<const Epetra_FECrsMatrix>
KeoFactory::
getKeo() const
{
    if ( !Ginla::Helpers::locaParameterVectorsEqual( keoBuildParameters_, mvp_->getParameters() ) )
    {
        this->fillKeo_( keo_, MATRIX_TYPE_REGULAR );
        keoBuildParameters_ = Teuchos::rcp( mvp_->getParameters()->clone() );
    }
    return keo_;
}
// =============================================================================
Teuchos::RCP<const Epetra_FECrsMatrix>
KeoFactory::
getKeoDMu() const
{
    if ( !Ginla::Helpers::locaParameterVectorsEqual( keoDMuBuildParameters_, mvp_->getParameters() ) )
    {
        this->fillKeo_( keoDMu_, MATRIX_TYPE_DMU );
        keoDMuBuildParameters_ = Teuchos::rcp( mvp_->getParameters()->clone() );
    }
    return keoDMu_;
}
// =============================================================================
Teuchos::RCP<const Epetra_FECrsMatrix>
KeoFactory::
getKeoDTheta() const
{
    if ( !Ginla::Helpers::locaParameterVectorsEqual( keoDThetaBuildParameters_, mvp_->getParameters() ) )
    {
        this->fillKeo_( keoDTheta_, MATRIX_TYPE_DTHETA );
        keoDThetaBuildParameters_ = Teuchos::rcp( mvp_->getParameters()->clone() );
    }
    return keoDTheta_;
}
// =============================================================================
const Teuchos::RCP<Epetra_FECrsGraph>
KeoFactory::
buildKeoGraph_() const
{
  // Which row/column map to use for the matrix?
  // The two possibilites are the non-overlapping map fetched from
  // the ownedNodes map, and the overlapping one from the
  // overlapNodes.
  // Let's illustrate the implications with the example of the matrix
  //   [ 2 1   ]
  //   [ 1 2 1 ]
  //   [   1 2 ].
  // Suppose subdomain 1 consists of node 1, subdomain 2 of node 3,
  // and node 2 forms the boundary between them.
  // For two processes, if process 1 owns nodes 1 and 2, the matrix
  // will be split as
  //   [ 2 1   ]   [       ]
  //   [ 1 2 1 ] + [       ]
  //   [       ]   [   1 2 ].
  // The vectors always need to have a unique map (otherwise, norms
  // cannot be computed by Epetra), so let's assume they have the
  // map ( [1,2], [3] ).
  // The communucation for a matrix-vector multiplication Ax=y
  // needs to be:
  //
  //   1. Communicate x(3) to process 1.
  //   2. Communicate x(2) to process 2.
  //   3. Compute.
  //
  // If the matrix is split up like
  //   [ 2 1   ]   [       ]
  //   [ 1 1   ] + [   1 1 ]
  //   [       ]   [   1 2 ]
  // (like the overlap map suggests), then any Ax=y comes down to:
  //
  //   1. Communicate x(2) to process 2.
  //   2. Compute.
  //   3. Communicate (part of) y(2) to process 1.
  //
  // In the general case, assuming that the number of nodes adjacent
  // to a boundary (on one side) are approximately the number of
  // nodes on that boundary, there is not much difference in
  // communication between the patterns.
  // What does differ, though, is the workload on the processes
  // during the computation phase: Process 1 that owns the whole
  // boundary, has to compute more than process 2.
  // Notice, however, that the total number of computations is
  // lower in scenario 1 (7 vs. 8 FLOPs); the same is true for
  // storage.
  // Hence, it comes down to the question whether or not the
  // mesh generator provided a fair share of the boundary nodes.
  // If yes, then scenario 1 will yield approximately even
  // computation times; if not, then scenario 2 will guarantee
  // equal computation times at the cost of higher total
  // storage and computation needs.
  //
  // Remark:
  // This matrix will later be fed into ML. ML has certain restrictions as to
  // what maps can be used. One of those is that RowMatrixRowMap() and
  // OperatorRangeMap must be the same, and, if the matrix is square,
  // OperatorRangeMap and OperatorDomainMap must coincide too.
  //
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  const Epetra_Map & noMap = *mesh_->getComplexNonOverlapMap();
  Teuchos::RCP<Epetra_FECrsGraph> keoGraph
      = Teuchos::rcp( new Epetra_FECrsGraph( Copy, noMap, 0 ) );

  if ( mesh_->supportsEdges() )
      this->buildKeoGraphEdges_( keoGraph );
  else
      this->buildKeoGraphCellEdges_( keoGraph );

  // Make sure that domain and range map are non-overlapping (to make sure that
  // states psi can compute norms) and equal (to make sure that the matrix works
  // with ML).
  TEUCHOS_ASSERT_EQUALITY( 0, keoGraph->GlobalAssemble(noMap,noMap) );

  return keoGraph;
}
// =============================================================================
void
KeoFactory::
fillKeo_( const Teuchos::RCP<Epetra_FECrsMatrix> & keoMatrix,
          const EMatrixType matrixType
        ) const
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*keoFillTime_);
#endif
  // Zero-out the matrix
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->PutScalar( 0.0 ) );

#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !thickness_.is_null() );
  TEUCHOS_ASSERT( !mvp_.is_null() );
#endif
  if ( mesh_->supportsEdges() )
  {
      this->fillKeoEdges_( keoMatrix,
                           matrixType,
                           mesh_->getEdgeCoefficients()
                         );
  }
  else
  {
      this->fillKeoCellEdges_( keoMatrix,
                               matrixType,
                               mesh_->getEdgeCoefficientsFallback()
                             );
  }

  // calls FillComplete by default
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->GlobalAssemble() );
  return;
}
// =============================================================================
void
KeoFactory::
buildKeoGraphEdges_( const Teuchos::RCP<Epetra_FECrsGraph> & keoGraph ) const
{
  // Get owned edges. This might throw as edges are not available for shells
  // (aka surfaces).
  std::vector<stk::mesh::Entity*> edges = mesh_->getOverlapEdges();

  if ( !globalIndexCacheUpToDate_ )
      this->buildGlobalIndexCache_( edges );

  // Loop over all edges and put entries wherever two nodes are connected.
  for ( unsigned int k=0; k<edges.size(); k++ )
      TEUCHOS_ASSERT_EQUALITY( 0, keoGraph->InsertGlobalIndices( 4, globalIndexCache_[k].Values(),
                                                                 4, globalIndexCache_[k].Values()
                                                               )
                             );

  return;
}
// =============================================================================
void
KeoFactory::
fillKeoEdges_( const Teuchos::RCP<Epetra_FECrsMatrix> & keoMatrix,
               const EMatrixType matrixType,
               const Teuchos::ArrayRCP<const double> & edgeCoefficients
             ) const
{
  // Get owned edges. This might throw as edges are not available for shells
  // (aka surfaces).
  const std::vector<stk::mesh::Entity*> edges = mesh_->getOverlapEdges();

  if ( !globalIndexCacheUpToDate_ )
      this->buildGlobalIndexCache_( edges );

  if ( !alphaCacheUpToDate_ )
      this->buildAlphaCache_( edges, edgeCoefficients );

  // Loop over all edges.
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
      // ---------------------------------------------------------------
      // Compute the integral
      //
      //    I = \int_{x0}^{xj} (xj-x0).A(x) / |xj-x0| dx
      //
      // numerically by the midpoint rule, i.e.,
      //
      //    I ~ |xj-x0| * (xj-x0) . A( 0.5*(xj+x0) ) / |xj-x0|.
      //
      // -------------------------------------------------------------------
      // Project vector field onto the edge.
      // Instead of first computing the projection over the normalized edge
      // and then multiply it with the edge length, don't normalize the
      // edge vector.
      double aInt = mvp_->getAEdgeMidpointProjection( k );

      double c, s, d;
      double sinAInt, cosAInt;
      //sinAInt = sin(aInt); cosAInt = cos(aInt);
      sincos( aInt, &sinAInt, &cosAInt );
      // For a slight speedup, move this switch statement out of the loop.
      switch ( matrixType )
      {
          case MATRIX_TYPE_REGULAR: // no derivative
          {
              c = alphaCache_[k] * cosAInt;
              s = alphaCache_[k] * sinAInt;
              d = - alphaCache_[k];
              break;
          }
          case MATRIX_TYPE_DMU: // dK/dmu
          {
              double dAdMuInt = mvp_->getdAdMuEdgeMidpointProjection( k );
              c = - alphaCache_[k] * dAdMuInt * sinAInt;
              s =   alphaCache_[k] * dAdMuInt * cosAInt;
              d = 0.0;
              break;
          }
          case MATRIX_TYPE_DTHETA: // dK/dtheta
          {
              double dAdThetaInt = mvp_->getdAdThetaEdgeMidpointProjection( k );
              c = - alphaCache_[k] * dAdThetaInt * sinAInt;
              s =   alphaCache_[k] * dAdThetaInt * cosAInt;
              d = 0.0;
              break;
          }
          default:
              TEST_FOR_EXCEPT_MSG( true,
                                    "Illegal matrix type \"" << matrixType << "\"."
                                  );
      }
      // We'd like to insert the 2x2 matrix
      //
      //     [   alpha                   , - alpha * exp( -IM * aInt ) ]
      //     [ - alpha * exp( IM * aInt ),   alpha                       ]
      //
      // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
      // that shares and edge.
      // Do that now, just blockwise for real and imaginary part.
      Epetra_SerialDenseMatrix A(4,4);
      A(0,0) = d;
      A(0,1) = 0.0;
      A(0,2) = c;
      A(0,3) = s;
      A(1,0) = 0.0;
      A(1,1) = d;
      A(1,2) = -s;
      A(1,3) = c;
      A(2,0) = c;
      A(2,1) = -s;
      A(2,2) = d;
      A(2,3) = 0.0;
      A(3,0) = s;
      A(3,1) = c;
      A(3,2) = 0.0;
      A(3,3) = d;
      TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues( globalIndexCache_[k], A ) );
      // -------------------------------------------------------------------
  }

  return;
}
// =============================================================================
void
KeoFactory::
buildGlobalIndexCache_( const std::vector<stk::mesh::Entity*> & edges ) const
{
  globalIndexCache_ = Teuchos::ArrayRCP<Epetra_IntSerialDenseVector>(edges.size());

  Teuchos::Tuple<int,2> gid;
  const stk::mesh::EntityRank nodeRank = mesh_->getMetaData()->node_rank();
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
      stk::mesh::PairIterRelation endPoints;
      endPoints = edges[k]->relations( nodeRank );

      gid[0] = (*endPoints[0].entity()).identifier() - 1;
      gid[1] = (*endPoints[1].entity()).identifier() - 1;

      globalIndexCache_[k] = Epetra_IntSerialDenseVector(4);
      globalIndexCache_[k][0] = 2*gid[0];
      globalIndexCache_[k][1] = 2*gid[0]+1;
      globalIndexCache_[k][2] = 2*gid[1];
      globalIndexCache_[k][3] = 2*gid[1]+1;
  }

  globalIndexCacheUpToDate_ = true;

  return;
}
// =============================================================================
void
KeoFactory::
buildAlphaCache_( const std::vector<stk::mesh::Entity*> & edges,
                  const Teuchos::ArrayRCP<const double> & edgeCoefficients
                ) const
{
  alphaCache_ = Teuchos::ArrayRCP<double>( edges.size() );

  Teuchos::Tuple<int,2> gid;
  const stk::mesh::EntityRank nodeRank = mesh_->getMetaData()->node_rank();
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
      stk::mesh::PairIterRelation endPoints;
      endPoints = edges[k]->relations( nodeRank );

      gid[0] = (*endPoints[0].entity()).identifier() - 1;
      gid[1] = (*endPoints[1].entity()).identifier() - 1;

      int tlid0 = thickness_->Map().LID( gid[0] );
      int tlid1 = thickness_->Map().LID( gid[1] );
#ifdef _DEBUG_
      TEST_FOR_EXCEPT_MSG( tlid0 < 0,
                           "The global index " << gid[0]
                           << " does not seem to be present on this node." );
      TEST_FOR_EXCEPT_MSG( tlid1 < 0,
                           "The global index " << gid[1]
                           << " does not seem to be present on this node." );
#endif
      double thickness = 0.5 * ( (*thickness_)[tlid0] + (*thickness_)[tlid1] );

      alphaCache_[k] = edgeCoefficients[k] * thickness;
  }

  alphaCacheUpToDate_ = true;

  return;
}
// =============================================================================
void
KeoFactory::
buildKeoGraphCellEdges_( const Teuchos::RCP<Epetra_FECrsGraph> & keoGraph ) const
{
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  if (!globalIndexFallbackCacheUpToDate_)
      this->buildGlobalIndexFallbackCache_( cells );

  // The index sets per edge per cell have been created in
  // buildGlobalIndexFallbackCache_().
  // Now, just insert them all; the order doesn't matter.
  for ( unsigned int k=0; k<cells.size(); k++ )
      for ( unsigned int l=0; l<globalIndexFallbackCache_[k].size(); l++ )
          TEUCHOS_ASSERT_EQUALITY( 0, keoGraph->InsertGlobalIndices( 4, globalIndexFallbackCache_[k][l].Values(),
                                                                     4, globalIndexFallbackCache_[k][l].Values()
                                                                   )
                                 );

  return;
}
// =============================================================================
void
KeoFactory::
fillKeoCellEdges_( const Teuchos::RCP<Epetra_FECrsMatrix> & keoMatrix,
                   const EMatrixType matrixType,
                   const Teuchos::ArrayRCP<const DoubleVector> & edgeCoefficientsFallback
                 ) const
{
  // get owned cells
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  if (!globalIndexFallbackCacheUpToDate_)
      this->buildGlobalIndexFallbackCache_( cells );

  if (!alphaFallbackCacheUpToDate_)
      this->buildAlphaFallbackCache_( cells, edgeCoefficientsFallback );

  // Loop over all edges.
  // To this end, loop over all cells and the edges within the cell.
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // Assume that the number of edges is correctly given by
      // globalIndexFallbackCache_.
      unsigned int numLocalEdges = globalIndexFallbackCache_[k].size();
      for ( unsigned int edgeIndex=0; edgeIndex < numLocalEdges; edgeIndex++ )
      {
          // -------------------------------------------------------------------
          // Compute the integral
          //
          //    I = \int_{x0}^{xj} (xj-x0).A(x) / |xj-x0| dx
          //
          // numerically by the midpoint rule, i.e.,
          //
          //    I ~ |xj-x0| * (xj-x0) . A( 0.5*(xj+x0) ) / |xj-x0|.
          //
          // -------------------------------------------------------------------
          // Project vector field onto the edge.
          // Instead of first computing the projection over the normalized edge
          // and then multiply it with the edge length, don't normalize the
          // edge vector.
          double aInt = mvp_->getAEdgeMidpointProjectionFallback( k, edgeIndex );
          double alpha = alphaFallbackCache_[k][edgeIndex];
          double c, s, d;
          double sinAInt, cosAInt;
          //sinAInt = sin(aInt);
          //cosAInt = cos(aInt);
          sincos( aInt, &sinAInt, &cosAInt );
          switch ( matrixType )
          {
              case MATRIX_TYPE_REGULAR: // no derivative
              {
                  c = alpha * cosAInt;
                  s = alpha * sinAInt;
                  d = - alpha;
                  break;
              }
              case MATRIX_TYPE_DMU: // dK/dmu
              {
                  double dAdMuInt = mvp_->getdAdMuEdgeMidpointProjectionFallback( k, edgeIndex );
                  c = - alpha * dAdMuInt * sinAInt;
                  s =   alpha * dAdMuInt * cosAInt;
                  d = 0.0;
                  break;
              }
              case MATRIX_TYPE_DTHETA: // dK/dtheta
              {
                  double dAdThetaInt = mvp_->getdAdThetaEdgeMidpointProjectionFallback( k, edgeIndex );
                  c = - alpha * dAdThetaInt * sinAInt;
                  s =   alpha * dAdThetaInt * cosAInt;
                  d = 0.0;
                  break;
              }
              default:
                  TEST_FOR_EXCEPT_MSG( true,
                                       "Illegal matrix type \"" << matrixType << "\"."
                                     );
          }

          // We'd like to insert the 2x2 matrix
          //
          //     [   alpha                   , - alpha * exp( -IM * aInt ) ]
          //     [ - alpha * exp( IM * aInt ),   alpha                       ]
          //
          // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
          // that shares and edge.
          // Do that now, just blockwise for real and imaginary part.
          Epetra_SerialDenseMatrix A(4,4);
          A(0,0) = d;
          A(0,1) = 0.0;
          A(0,2) = c;
          A(0,3) = s;
          A(1,0) = 0.0;
          A(1,1) = d;
          A(1,2) = -s;
          A(1,3) = c;
          A(2,0) = c;
          A(2,1) = -s;
          A(2,2) = d;
          A(2,3) = 0.0;
          A(3,0) = s;
          A(3,1) = c;
          A(3,2) = 0.0;
          A(3,3) = d;
          TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues( globalIndexFallbackCache_[k][edgeIndex], A ) );
          // -------------------------------------------------------------------
      }
  }

  return;
}
// =============================================================================
void
KeoFactory::
buildGlobalIndexFallbackCache_( const std::vector<stk::mesh::Entity*> & cells ) const
{
  globalIndexFallbackCache_ =
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<Epetra_IntSerialDenseVector> >(cells.size());

  int gid0, gid1;
  Epetra_IntSerialDenseVector gIndices(4);
  const stk::mesh::EntityRank nodeRank = mesh_->getMetaData()->node_rank();
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation nodes = cells[k]->relations( nodeRank );
      unsigned int numLocalNodes = nodes.size();
      unsigned int numLocalEdges = numLocalNodes*(numLocalNodes-1) / 2;
      unsigned int edgeIndex = 0;

      globalIndexFallbackCache_[k] =
          Teuchos::ArrayRCP<Epetra_IntSerialDenseVector>(numLocalEdges);

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
      {
          gid0 = (*nodes[e0].entity()).identifier() - 1;
          gIndices[0] = 2*gid0;
          gIndices[1] = 2*gid0+1;
          for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
          {
              gid1 = (*nodes[e1].entity()).identifier() - 1;
              gIndices[2] = 2*gid1;
              gIndices[3] = 2*gid1+1;
              globalIndexFallbackCache_[k][edgeIndex] = gIndices;
              edgeIndex++;
          }
      }
  }

  globalIndexFallbackCacheUpToDate_ = true;

  return;
}
// =============================================================================
void
KeoFactory::
buildAlphaFallbackCache_( const std::vector<stk::mesh::Entity*> & cells,
                          const Teuchos::ArrayRCP<const DoubleVector> & edgeCoefficientsFallback
                        ) const
{
  unsigned int numCells = cells.size();
  alphaFallbackCache_ = Teuchos::ArrayRCP<DoubleVector>( numCells );

  Teuchos::Tuple<int,2> gid;
  const stk::mesh::EntityRank nodeRank = mesh_->getMetaData()->node_rank();
  for ( unsigned int k=0; k<numCells; k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation nodes = cells[k]->relations( nodeRank );

      unsigned int numLocalNodes = nodes.size();
      unsigned int numLocalEdges = numLocalNodes * (numLocalNodes-1) / 2;
      alphaFallbackCache_[k] = DoubleVector(numLocalEdges);
      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      unsigned int edgeIndex = 0;
      int gid0, gid1;
      for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
      {
          gid0 = (*nodes[e0].entity()).identifier() - 1;
          int tlid0 = thickness_->Map().LID( gid0 );
#ifdef _DEBUG_
          TEST_FOR_EXCEPT_MSG( tlid0 < 0,
                               "The global index " << gid0
                               << " does not seem to be present on this node." );
#endif
          for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
          {
              gid1 = (*nodes[e1].entity()).identifier() - 1;
              int tlid1 = thickness_->Map().LID( gid1 );
#ifdef _DEBUG_
              TEST_FOR_EXCEPT_MSG( tlid1 < 0,
                                   "The global index " << gid1
                                   << " does not seem to be present on this node." );
#endif
              double thickness = 0.5 * ( (*thickness_)[tlid0] + (*thickness_)[tlid1] );

              // Multiply by the thickness value of the midpoint. As this is not
              // available, take the mean between the values at the nodes.
              alphaFallbackCache_[k][edgeIndex] = edgeCoefficientsFallback[k][edgeIndex]
                                                * thickness;

              edgeIndex++;
              // -------------------------------------------------------------------
          }
      }
  }

  alphaFallbackCacheUpToDate_ = true;

  return;
}
// =============================================================================
} // namespace Ginla
