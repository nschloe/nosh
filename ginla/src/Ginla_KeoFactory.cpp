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
// includes
#include "Ginla_KeoFactory.hpp"
#include "Ginla_StkMesh.hpp"
#include "Ginla_Helpers.hpp"

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
            const Teuchos::RCP<Ginla::MagneticVectorPotential> & mvp
          ):
#ifdef GINLA_TEUCHOS_TIME_MONITOR
        buildKeoTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo") ),
        buildKeoGraphTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeoGraph") ),
        buildKeoTime1_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo01") ),
        buildKeoTime2_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo02") ),
        buildKeoTime3_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo03") ),
        buildKeoTime4_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo04") ),
        buildKeoTime5_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo05") ),
        buildKeoTime6_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo06") ),
        buildKeoTime7_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo07") ),
        buildKeoTime8_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo08") ),
        buildKeoTime9_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo09") ),
        buildKeoTime10_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo10") ),
        buildKeoTime11_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo11") ),
        buildKeoTime12_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo12") ),
#endif
        mesh_ ( mesh ),
        thickness_( thickness ),
        mvp_( mvp ),
        keoGraph_( this->buildKeoGraph_() ), // build the graph immediately
        keo_( Teuchos::rcp( new Epetra_CrsMatrix( Copy, *keoGraph_ ) ) ),
        keoBuildParameters_( Teuchos::null ),
        keoDMu_( Teuchos::rcp( new Epetra_CrsMatrix( Copy, *keoGraph_ ) ) ),
        keoDMuBuildParameters_( Teuchos::null ),
        localColIndexCache_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> >() ),
        localRowIndexCache_( Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> >() ),
        localIndexCacheUpToDate_( false ),
        alphaCache_( Teuchos::ArrayRCP<double>() ),
        alphaCacheUpToDate_( false )
{
}
// =============================================================================
KeoFactory::
~KeoFactory()
{
}
// =============================================================================
const Teuchos::RCP<const Ginla::StkMesh>
KeoFactory::
getMesh() const
{
    return mesh_;
}
// =============================================================================
void
KeoFactory::
updateParameters( const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams
                ) const
{
  // set the parameters
  TEUCHOS_ASSERT( !mvpParams.is_null() );
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
void
KeoFactory::
fillKeo_( const Teuchos::RCP<Epetra_CrsMatrix> & keoMatrix,
          const EMatrixType matrixType
        ) const
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime_);
#endif
  // Zero-out the matrix
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->PutScalar( 0.0 ) );

  Teuchos::ArrayRCP<DoubleVector> edgeCoefficients;
  Teuchos::ArrayRCP<DoubleVector> edgeCoefficientsFallback;
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime1_);
#endif
  TEUCHOS_ASSERT( !mesh_.is_null() );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.

  TEUCHOS_ASSERT( !thickness_.is_null() );
  TEUCHOS_ASSERT( !mvp_.is_null() );

  try
  {
      this->fillKeoEdges_( keoMatrix,
                           matrixType,
                           mesh_->getEdgeCoefficients()
                         );
  }
  catch( ... )
  {
      this->fillKeoCellEdges_( keoMatrix,
                               matrixType,
                               mesh_->getEdgeCoefficientsFallback()
                             );
  }

}

{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime12_);
#endif
  // calls FillComplete by default
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->FillComplete() );
}
  return;
}
// =============================================================================
void
KeoFactory::
fillKeoEdges_( const Teuchos::RCP<Epetra_CrsMatrix> & keoMatrix,
               const EMatrixType matrixType,
               const Teuchos::ArrayRCP<const double> & edgeCoefficients
             ) const
{
  // get owned edges
  const std::vector<stk::mesh::Entity*> edges = mesh_->getOverlapEdges();

  if ( !localIndexCacheUpToDate_ )
      this->buildLocalIndexCache_( edges );

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
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime6_);
#endif
      double v[4];
      // sum it all in!
      v[0] = d;
      v[1] = 0.0;
      v[2] = c;
      v[3] = s;
      TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoMyValues( localRowIndexCache_[k][0], 4, v, localColIndexCache_[k].getRawPtr() ) );
      v[0] = 0.0;
      v[1] = d;
      v[2] = -s;
      v[3] = c;
      TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoMyValues( localRowIndexCache_[k][1], 4, v, localColIndexCache_[k].getRawPtr() ) );
      v[0] = c;
      v[1] = -s;
      v[2] = d;
      v[3] = 0.0;
      TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoMyValues( localRowIndexCache_[k][2], 4, v, localColIndexCache_[k].getRawPtr() ) );
      v[0] = s;
      v[1] = c;
      v[2] = 0.0;
      v[3] = d;
      TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoMyValues( localRowIndexCache_[k][3], 4, v, localColIndexCache_[k].getRawPtr() ) );
}
      // -------------------------------------------------------------------
  }

  return;
}
// =============================================================================
void
KeoFactory::
buildLocalIndexCache_( const std::vector<stk::mesh::Entity*> & edges ) const
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime2_);
#endif

  localColIndexCache_ = Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> >(edges.size());
  localRowIndexCache_ = Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> >(edges.size());

  Teuchos::Tuple<int,2> gid;
  Teuchos::Tuple<int,2> lid;
  int gIndices[4];
  int localColIndices[4];
  int localRowIndices[4];
  const stk::mesh::EntityRank nodeRank = mesh_->getMetaData()->node_rank();
  for ( unsigned int k=0; k<edges.size(); k++ )
  {
      stk::mesh::PairIterRelation endPoints;
      endPoints = edges[k]->relations( nodeRank );

      gid[0] = (*endPoints[0].entity()).identifier() - 1;
      gid[1] = (*endPoints[1].entity()).identifier() - 1;

      gIndices[0] = 2*gid[0];
      gIndices[1] = 2*gid[0]+1;
      gIndices[2] = 2*gid[1];
      gIndices[3] = 2*gid[1]+1;

      localColIndexCache_[k] = Teuchos::ArrayRCP<int>(4);
      localColIndexCache_[k][0] = keoGraph_->ColMap().LID( gIndices[0] );
//       TEST_FOR_EXCEPT_MSG( localColIndices[0] < 0,
//                           "The global index " << gIndices[0]
//                           << " does not seem to be present on this node." );
      localColIndexCache_[k][1] = keoGraph_->ColMap().LID( gIndices[1] );
//       TEST_FOR_EXCEPT_MSG( localColIndices[0] < 0,
//                           "The global index " << gIndices[1]
//                           << " does not seem to be present on this node." );
      localColIndexCache_[k][2] = keoGraph_->ColMap().LID( gIndices[2] );
//       TEST_FOR_EXCEPT_MSG( localColIndices[2] < 0,
//                           "The global index " << gIndices[2]
//                           << " does not seem to be present on this node." );
      localColIndexCache_[k][3] = keoGraph_->ColMap().LID( gIndices[3] );
//       TEST_FOR_EXCEPT_MSG( localColIndices[3] < 0,
//                           "The global index " << gIndices[3]
//                           << " does not seem to be present on this node." );

      localRowIndexCache_[k] = Teuchos::ArrayRCP<int>(4);
      localRowIndexCache_[k][0] = keoGraph_->RowMap().LID( gIndices[0] );
//       TEST_FOR_EXCEPT_MSG( localRowIndices[0] < 0,
//                           "The global index " << gIndices[0]
//                           << " does not seem to be present on this node." );
      localRowIndexCache_[k][1] = keoGraph_->RowMap().LID( gIndices[1] );
//       TEST_FOR_EXCEPT_MSG( localRowIndices[0] < 0,
//                           "The global index " << gIndices[1]
//                           << " does not seem to be present on this node." );
      localRowIndexCache_[k][2] = keoGraph_->RowMap().LID( gIndices[2] );
//       TEST_FOR_EXCEPT_MSG( localRowIndices[2] < 0,
//                           "The global index " << gIndices[2]
//                           << " does not seem to be present on this node." );
      localRowIndexCache_[k][3] = keoGraph_->RowMap().LID( gIndices[3] );
//       TEST_FOR_EXCEPT_MSG( localRowIndices[3] < 0,
//                           "The global index " << gIndices[3]
//                           << " does not seem to be present on this node." );

  }

  localIndexCacheUpToDate_ = true;

  return;
}
// =============================================================================
void
KeoFactory::
buildAlphaCache_( const std::vector<stk::mesh::Entity*> & edges,
                  const Teuchos::ArrayRCP<const double> & edgeCoefficients
                ) const
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime3_);
#endif

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
      TEST_FOR_EXCEPT_MSG( tlid0 < 0,
                           "The global index " << gid[0]
                           << " does not seem to be present on this node." );
      int tlid1 = thickness_->Map().LID( gid[1] );
      TEST_FOR_EXCEPT_MSG( tlid1 < 0,
                           "The global index " << gid[1]
                           << " does not seem to be present on this node." );
      double thickness = 0.5 * ( (*thickness_)[tlid0] + (*thickness_)[tlid1] );

      alphaCache_[k] = edgeCoefficients[k] * thickness;
  }

  alphaCacheUpToDate_ = true;

  return;
}
// =============================================================================
void
KeoFactory::
fillKeoCellEdges_( const Teuchos::RCP<Epetra_CrsMatrix> & keoMatrix,
                   const EMatrixType matrixType,
                   const Teuchos::ArrayRCP<const DoubleVector> & edgeCoefficientsFallback
                 ) const
{
  // get owned cells
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // Loop over all edges.
  // To this end, loop over all cells and the edges within the cell.
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      unsigned int numLocalNodes;
      // extract the nodal coordinates
      stk::mesh::PairIterRelation rel;
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime2_);
#endif
      // get the nodes local to the cell
      rel = cells[k]->relations( mesh_->getMetaData()->node_rank() );

      numLocalNodes = rel.size();
}
      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      unsigned int edgeIndex = 0;
      Teuchos::Tuple<int,2> gid;
      Teuchos::Tuple<int,2> lid;
      int gIndices[4];
      int lIndices[4];
      for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
      {
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime3_);
#endif
          gid[0] = (*rel[e0].entity()).identifier() - 1;
          gIndices[0] = 2*gid[0];
          gIndices[1] = 2*gid[0]+1;
          lIndices[0] = keoGraph_->ColMap().LID( gIndices[0] );
//          TEST_FOR_EXCEPT_MSG( lIndices[0] < 0,
//                               "The global index " << gIndices[0]
//                               << " does not seem to be present on this node." );
          lIndices[1] = keoGraph_->ColMap().LID( gIndices[1] );
//          TEST_FOR_EXCEPT_MSG( lIndices[0] < 0,
//                               "The global index " << gIndices[1]
//                               << " does not seem to be present on this node." );
}
          for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
          {
double aInt;
double alpha;
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime4_);
#endif
              gid[1] = (*rel[e1].entity()).identifier() - 1;
              gIndices[2] = 2*gid[1];
              gIndices[3] = 2*gid[1]+1;
              lIndices[2] = keoGraph_->ColMap().LID( gIndices[2] );
//              TEST_FOR_EXCEPT_MSG( lIndices[2] < 0,
//                                   "The global index " << gIndices[2]
//                                   << " does not seem to be present on this node." );
              lIndices[3] = keoGraph_->ColMap().LID( gIndices[3] );
//              TEST_FOR_EXCEPT_MSG( lIndices[3] < 0,
//                                   "The global index " << gIndices[3]
//                                   << " does not seem to be present on this node." );
}
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime5_);
#endif
              // edge weights
              alpha = edgeCoefficientsFallback[k][edgeIndex];
}
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime6_);
#endif
              // Multiply by the thickness value of the midpoint. As this is not available,
              // take the mean between the values at the nodes.
              int tlid0 = thickness_->Map().LID( gid[0] );
              TEST_FOR_EXCEPT_MSG( tlid0 < 0,
                                   "The global index " << gid[0]
                                   << " does not seem to be present on this node." );
              int tlid1 = thickness_->Map().LID( gid[1] );
              TEST_FOR_EXCEPT_MSG( tlid1 < 0,
                                   "The global index " << gid[1]
                                   << " does not seem to be present on this node." );
              double thickness = 0.5 * ( (*thickness_)[tlid0] + (*thickness_)[tlid1] );

              alpha *= thickness;
}
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
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime7_);
#endif
              // Filling the MVP cache takes another 50% of the time in this whole function.
// TODO move this in the outer loop?
              aInt = mvp_->getAEdgeMidpointProjectionFallback( k, edgeIndex );
}
              double c, s, d;
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime8_);
#endif
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
                  default:
                      TEST_FOR_EXCEPT_MSG( true,
                                            "Illegal matrix type \"" << matrixType << "\"."
                                          );
              }

}
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime9_);
#endif

              // We'd like to insert the 2x2 matrix
              //
              //     [   alpha                   , - alpha * exp( -IM * aInt ) ]
              //     [ - alpha * exp( IM * aInt ),   alpha                       ]
              //
              // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
              // that shares and edge.
              // Do that now, just blockwise for real and imaginary part.
              //indices[2] = 2*lid[1];
              //indices[3] = 2*lid[1]+1;

}
              double v[4];
              int ii[3];
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime10_);
#endif
              // sum it all in!
//              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues ( indices, values ) );
              v[0] = d;
              v[1] = 0.0;
              v[2] = c;
              v[3] = s;
/*              ii[0] = lIndices[0];
              ii[1] = lIndices[2];
              ii[2] = lIndices[3];*/
}
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime11_);
#endif
              int ddd = keoMatrix->RowMap().LID( gIndices[0] );
              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoMyValues( ddd, 4, v, lIndices ) );
//              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues( gIndices[0], 4, v, gIndices ) );
}
              v[0] = 0.0;
              v[1] = d;
              v[2] = -s;
              v[3] = c;
              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues( gIndices[1], 4, v, gIndices ) );
              v[0] = c;
              v[1] = -s;
              v[2] = d;
              v[3] = 0.0;
              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues( gIndices[2], 4, v, gIndices ) );
              v[0] = s;
              v[1] = c;
              v[2] = 0.0;
              v[3] = d;
              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues( gIndices[3], 4, v, gIndices ) );

              edgeIndex++;
              // -------------------------------------------------------------------
          }
      }
  }

  return;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
KeoFactory::
buildKeo_( const EMatrixType matrixType ) const
{
  TEUCHOS_ASSERT( !keoGraph_.is_null() );

  // create an zeroed-out matrix
  Teuchos::RCP<Epetra_CrsMatrix> keoMatrix =
      Teuchos::rcp( new Epetra_CrsMatrix( Copy, *keoGraph_ ) );

  // Fill the matrix.
  this->fillKeo_( keoMatrix, matrixType );

  return keoMatrix;
}
// =============================================================================
Teuchos::RCP<const Epetra_CrsMatrix>
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
Teuchos::RCP<const Epetra_CrsMatrix>
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
Teuchos::RCP<const Epetra_CrsGraph>
KeoFactory::
getKeoGraph() const
{
  TEUCHOS_ASSERT( !keoGraph_.is_null() );
  return keoGraph_;
}
// =============================================================================
const Teuchos::RCP<Epetra_CrsGraph>
KeoFactory::
buildKeoGraph_() const
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoGraphTime_);
#endif
  // This method used to be implemented with Epetra_FECrsGraph,
  // but it turns out that there's no gain in runtime performance
  // compare to a regular Epetra_CrsGraph implementation.

  TEUCHOS_ASSERT( !mesh_.is_null() );

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
  Teuchos::RCP<Epetra_CrsGraph> keoGraph
      = Teuchos::rcp( new Epetra_CrsGraph( Copy, *mesh_->getComplexOverlapMap(), 0 ) );
//  Teuchos::RCP<Epetra_CrsGraph> keoGraph
//      = Teuchos::rcp( new Epetra_CrsGraph( Copy, *mesh_->getComplexNonOverlapMap(), 0 ) );

  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // Loop over all edges and put entries wherever two nodes are connected.
  Teuchos::Tuple<int,2> gid;
//  Teuchos::Tuple<int,2> lid;
  int indices[4];
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation rel = cells[k]->relations( mesh_->getMetaData()->node_rank() );

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      for ( unsigned int e0=0; e0<rel.size(); e0++ )
      {
          gid[0] = (*rel[e0].entity()).identifier() - 1;
//          lid[0] = mesh_->getNodesOverlapMap()->LID( gid[0] );
//          TEST_FOR_EXCEPT_MSG( lid[0] < 0,
//                               "The global index " << gid[0]
//                               << " does not seem to be present on this node ("
//                               << mesh_->getComm().MyPID()
//                               << ")." );
          indices[0] = 2*gid[0];
          indices[1] = 2*gid[0]+1;
          for ( unsigned int e1=e0+1; e1<rel.size(); e1++ )
          {
              gid[1] = (*rel[e1].entity()).identifier() - 1;
//              lid[1] = mesh_->getNodesOverlapMap()->LID( gid[1] );
//              TEST_FOR_EXCEPT_MSG( lid[1] < 0,
//                                   "The global index " << gid[1]
//                                   << " does not seem to be present on this node ("
//                                   << mesh_->getComm().MyPID()
//                                   << ")." );
              indices[2] = 2*gid[1];
              indices[3] = 2*gid[1]+1;

              // We cannot use InsertMyIndices here as FillComplete() hasn't been called yet
              // and hence IndicesAreLocal()==false.
              TEUCHOS_ASSERT_EQUALITY( 0, keoGraph->InsertGlobalIndices( indices[0], 4, indices ) );
              TEUCHOS_ASSERT_EQUALITY( 0, keoGraph->InsertGlobalIndices( indices[1], 4, indices ) );
              TEUCHOS_ASSERT_EQUALITY( 0, keoGraph->InsertGlobalIndices( indices[2], 4, indices ) );
              TEUCHOS_ASSERT_EQUALITY( 0, keoGraph->InsertGlobalIndices( indices[3], 4, indices ) );
          }
      }
  }
  // Make sure that domain and range map are non-overlapping.
  const Epetra_Map & noMap = *mesh_->getComplexNonOverlapMap();
//  TEUCHOS_ASSERT_EQUALITY( 0, keoGraph.GlobalAssemble(noMap,noMap) );
  TEUCHOS_ASSERT_EQUALITY( 0, keoGraph->FillComplete(noMap,noMap) );

  return keoGraph;
}
// =============================================================================
} // namespace Ginla
