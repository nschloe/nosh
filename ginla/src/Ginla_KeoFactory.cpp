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

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>

#include <ml_epetra_preconditioner.h>

#include <Teuchos_TimeMonitor.hpp>

namespace Ginla {
// =============================================================================
KeoFactory::
KeoFactory( const Teuchos::RCP<Ginla::StkMesh>      & mesh,
            const Teuchos::RCP<const Epetra_Vector>            & thickness,
            const Teuchos::RCP<Ginla::MagneticVectorPotential> & mvp
          ):
        mesh_ ( mesh ),
        thickness_( thickness ),
        mvp_( mvp ),
        buildKeoTime_( Teuchos::TimeMonitor::getNewTimer("KeoFactory::buildKeo") )
{
}
// =============================================================================
KeoFactory::
~KeoFactory()
{
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
buildKeo( Epetra_FECrsMatrix & keoMatrix ) const
{
  Teuchos::TimeMonitor tm(*buildKeoTime_);

  // zero out the matrix
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix.PutScalar( 0.0 ) );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
  // get owned cells
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // This takes about 50% of the time in this whole function.
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeCoefficients = mesh_->getEdgeCoefficients();
  TEUCHOS_ASSERT( !edgeCoefficients.is_null() );

  // Loop over all edges.
  // To this end, loop over all cells and the edges within the cell.
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation rel = (*cells[k]).relations();

      unsigned int numLocalNodes = rel.size();
      // extract the nodal coordinates
      Teuchos::Array<Point> localNodes;

      localNodes = mesh_->getNodeCoordinates( rel );

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      unsigned int edgeIndex = 0;
      Teuchos::Tuple<int,2> gid;
      Teuchos::Tuple<int,2> lid;
      for ( unsigned int e0 = 0; e0 < numLocalNodes; e0++ )
      {
          const Point & node0 = localNodes[e0];
          gid[0] = (*rel[e0].entity()).identifier() - 1;
          lid[0] = thickness_->Map().LID( gid[0] );
          TEST_FOR_EXCEPT_MSG( lid[0] < 0,
                                       "The global index " << gid[0]
                                       << " does not seem to be present on this node." );
          for ( unsigned int e1 = e0+1; e1 < numLocalNodes; e1++ )
          {
              const Point & node1 = localNodes[e1];
              gid[1] = (*rel[e1].entity()).identifier() - 1;
              lid[1] = thickness_->Map().LID( gid[1] );
              TEST_FOR_EXCEPT_MSG( lid[1] < 0,
                                           "The global index " << gid[1]
                                           << " does not seem to be present on this node." );

              // coarea / edge ratio
              double alpha = edgeCoefficients[k][edgeIndex];

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
              double aInt = mvp_->getAEdgeMidpointProjection( k, edgeIndex );

              // We'd like to insert the 2x2 matrix
              //
              //     [   alpha                   , - alpha * exp( -IM * aInt ) ]
              //     [ - alpha * exp( IM * aInt ),   alpha                       ]
              //
              // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
              // that shares and edge.
              // Do that now, just blockwise for real and imaginary part.
              Epetra_IntSerialDenseVector indices(4);
              indices[0] = 2*gid[0];
              indices[1] = 2*gid[0]+1;
              indices[2] = 2*gid[1];
              indices[3] = 2*gid[1]+1;

              Epetra_SerialDenseMatrix values( 4, 4 );
              values(0,0) = - alpha;
              values(0,1) = 0.0;
              values(1,0) = 0.0;
              values(1,1) = - alpha;

              double alphaCosAInt = alpha * cos(aInt);
              double alphaSinAInt = alpha * sin(aInt);

              values(0,2) =   alphaCosAInt;
              values(0,3) =   alphaSinAInt;
              values(1,2) = - alphaSinAInt;
              values(1,3) =   alphaCosAInt;

              values(2,0) =   alphaCosAInt;
              values(2,1) = - alphaSinAInt;
              values(3,0) =   alphaSinAInt;
              values(3,1) =   alphaCosAInt;

              values(2,2) = - alpha;
              values(2,3) = 0.0;
              values(3,2) = 0.0;
              values(3,3) = - alpha;

              // sum it all in!
              TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix.SumIntoGlobalValues ( indices, values ) );

              edgeIndex++;
              // -------------------------------------------------------------------
          }
      }
  }

  // calls FillComplete by default
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix.GlobalAssemble() );

  return;
}
// =============================================================================
const Epetra_FECrsGraph
KeoFactory::
buildKeoGraph() const
{
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
  // and node 2 froms the boundary between them.
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
  const Epetra_Map & noMap = *mesh_->getComplexNonOverlapMap();
  Epetra_FECrsGraph keoGraph( Copy, noMap, 0 );
//  Epetra_FECrsGraph keoGraph( Copy, *mesh_->getComplexOverlapMap(), 0 );

  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  // Loop over all edges and put entries whereever two nodes are connected.
  Teuchos::Tuple<int,2> gid;
  for ( unsigned int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation rel = (*cells[k]).relations();

      // In a simplex, the edges are exactly the connection between each pair
      // of nodes. Hence, loop over pairs of nodes.
      for ( unsigned int e0=0; e0<rel.size(); e0++ )
      {
          gid[0] = (*rel[e0].entity()).identifier() - 1;
          for ( unsigned int e1=e0+1; e1<rel.size(); e1++ )
          {
              gid[1] = (*rel[e1].entity()).identifier() - 1;

              Epetra_IntSerialDenseVector indices(4);
              indices[0] = 2*gid[0];
              indices[1] = 2*gid[0]+1;
              indices[2] = 2*gid[1];
              indices[3] = 2*gid[1]+1;

              // sum it all in!
              TEUCHOS_ASSERT_EQUALITY( 0, keoGraph.InsertGlobalIndices( indices.Length(), indices.Values(), // rows
                                                                        indices.Length(), indices.Values()  // cols
                                                                      )
                                     );
          }
      }
  }

  // Make sure that domain and range map are non-overlapping.
  TEUCHOS_ASSERT_EQUALITY( 0, keoGraph.GlobalAssemble(noMap,noMap) );

  return keoGraph;
}
// =============================================================================
} // namespace Ginla
