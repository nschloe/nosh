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
        buildKeoTime1_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo1") ),
        buildKeoTime2_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo2") ),
        buildKeoTime3_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo3") ),
        buildKeoTime4_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo4") ),
        buildKeoTime5_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo5") ),
        buildKeoTime6_( Teuchos::TimeMonitor::getNewTimer("Ginla: KeoFactory::buildKeo5") ),
#endif
        mesh_ ( mesh ),
        thickness_( thickness ),
        mvp_( mvp ),
        keoGraph_( this->buildKeoGraph_() ), // build the graph immediately
        keo_( Teuchos::rcp( new Epetra_CrsMatrix( Copy, *keoGraph_ ) ) ),
        keoBuildParameters_( Teuchos::null ),
        keoDMu_( Teuchos::rcp( new Epetra_CrsMatrix( Copy, *keoGraph_ ) ) ),
        keoDMuBuildParameters_( Teuchos::null )
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
fillKeo_( const Teuchos::RCP<Epetra_CrsMatrix> keoMatrix,
          const EMatrixType matrixType
        ) const
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime_);
#endif
  // Zero-out the matrix
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->PutScalar( 0.0 ) );

  std::vector<stk::mesh::Entity*> cells;
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeCoefficients;
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime4_);
#endif
  TEUCHOS_ASSERT( !mesh_.is_null() );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
  // get owned cells
  cells = mesh_->getOwnedCells();

  // This takes about 50% of the time in this whole function.
  edgeCoefficients = mesh_->getEdgeCoefficients();
  TEUCHOS_ASSERT( !edgeCoefficients.is_null() );

  // Loop over all edges.
  // To this end, loop over all cells and the edges within the cell.
  TEUCHOS_ASSERT( !thickness_.is_null() );
  TEUCHOS_ASSERT( !mvp_.is_null() );
}

  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      unsigned int numLocalNodes;
      // extract the nodal coordinates
      Teuchos::Array<Point> localNodes;
      stk::mesh::PairIterRelation rel;
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime4_);
#endif
      // get the nodes local to the cell
      rel = (*cells[k]).relations();

      numLocalNodes = rel.size();

      localNodes = mesh_->getNodeCoordinates( rel );
}
      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      unsigned int edgeIndex = 0;
      Teuchos::Tuple<int,2> gid;
      Teuchos::Tuple<int,2> lid;
      int indices[4];
      for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
      {
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime5_);
#endif
          const Point & node0 = localNodes[e0];
          gid[0] = (*rel[e0].entity()).identifier() - 1;
          lid[0] = thickness_->Map().LID( gid[0] );
          indices[0] = 2*gid[0];
          indices[1] = 2*gid[0]+1;
          TEST_FOR_EXCEPT_MSG( lid[0] < 0,
                                       "The global index " << gid[0]
                                       << " does not seem to be present on this node." );
}
          for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
          {
double aInt;
double alpha;
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime1_);
#endif

              const Point & node1 = localNodes[e1];
              gid[1] = (*rel[e1].entity()).identifier() - 1;
              lid[1] = thickness_->Map().LID( gid[1] );
              TEST_FOR_EXCEPT_MSG( lid[1] < 0,
                                           "The global index " << gid[1]
                                           << " does not seem to be present on this node." );

              // edge weights
              alpha = edgeCoefficients[k][edgeIndex];

              // Multiply by the thickness value of the midpoint. As this is not available,
              // take the mean between the values at the nodes.
              double thickness = 0.5 * ( (*thickness_)[lid[0]] + (*thickness_)[lid[1]] );

              alpha *= thickness;

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

              // Filling the MVP cache takes another 50% of the time in this whole function.
              aInt = mvp_->getAEdgeMidpointProjection( k, edgeIndex );
}
              double c, s, d;
              switch ( matrixType )
              {
                  case MATRIX_TYPE_REGULAR: // no derivative
                  {
                      c = alpha * cos(aInt);
                      s = alpha * sin(aInt);
                      d = - alpha;
                      break;
                  }
                  case MATRIX_TYPE_DMU: // dK/dmu
                  {
                      double dAdMuInt = mvp_->getdAdMuEdgeMidpointProjection( k, edgeIndex );
                      c = - alpha * dAdMuInt * sin(aInt);
                      s =   alpha * dAdMuInt * cos(aInt);
                      d = 0.0;
                      break;
                  }
                  default:
                      TEST_FOR_EXCEPT_MSG( true,
                                           "Illegal matrix type \"" << matrixType << "\"."
                                         );
              }


{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime2_);
#endif

              // We'd like to insert the 2x2 matrix
              //
              //     [   alpha                   , - alpha * exp( -IM * aInt ) ]
              //     [ - alpha * exp( IM * aInt ),   alpha                       ]
              //
              // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
              // that shares and edge.
              // Do that now, just blockwise for real and imaginary part.
              indices[2] = 2*gid[1];
              indices[3] = 2*gid[1]+1;

}
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*buildKeoTime3_);
#endif
              // sum it all in!
//              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues ( indices, values ) );
              double v[4];
              v[0] = d;
              v[1] = 0.0;
              v[2] = c;
              v[3] = s;
              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues ( indices[0], 4, v, indices ) );
              v[0] = 0.0;
              v[1] = d;
              v[2] = -s;
              v[3] = c;
              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues ( indices[1], 4, v, indices ) );
              v[0] = c;
              v[1] = -s;
              v[2] = d;
              v[3] = 0.0;
              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues ( indices[2], 4, v, indices ) );
              v[0] = s;
              v[1] = c;
              v[2] = 0.0;
              v[3] = d;
              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->SumIntoGlobalValues ( indices[3], 4, v, indices ) );
}

              edgeIndex++;
              // -------------------------------------------------------------------
          }
      }
  }

  // calls FillComplete by default
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->FillComplete() );

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
  int indices[4];
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation rel = (*cells[k]).relations();

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      for ( unsigned int e0=0; e0<rel.size(); e0++ )
      {
          gid[0] = (*rel[e0].entity()).identifier() - 1;
          indices[0] = 2*gid[0];
          indices[1] = 2*gid[0]+1;
          for ( unsigned int e1=e0+1; e1<rel.size(); e1++ )
          {
              gid[1] = (*rel[e1].entity()).identifier() - 1;
              indices[2] = 2*gid[1];
              indices[3] = 2*gid[1]+1;

              // insert it all
//              TEUCHOS_ASSERT_EQUALITY( 0, keoGraph.InsertGlobalIndices( 4, indices, // rows
//                                                                        4, indices  // cols
//                                                                      )
//                                     );
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
