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
// includes
#include "Ginla_EpetraFVM_KeoFactory.h"

#include "Ginla_EpetraFVM_StkMesh.h"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosBlockCGSolMgr.hpp>

#include <ml_epetra_preconditioner.h>
// =============================================================================
Ginla::EpetraFVM::KeoFactory::
KeoFactory( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh>               & mesh,
            const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp
          ):
        mesh_ ( mesh ),
        mvp_( mvp )
{
}
// =============================================================================
Ginla::EpetraFVM::KeoFactory::
~KeoFactory()
{
}
// =============================================================================
void
Ginla::EpetraFVM::KeoFactory::
buildKeo( Epetra_FECrsMatrix                              & keoMatrix,
          const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams,
          const Teuchos::Tuple<double,3>                  & scaling
        ) const
{
  keoMatrix.PutScalar( 0.0 );

  // set the parameters
  TEUCHOS_ASSERT( !mvpParams.is_null() );
  mvp_->setParameters( *mvpParams );
  mesh_->scale( scaling );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the cells, create local load vector and mass matrix,
  // and insert them into the global matrix.
  // get owned cells
  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > coedgeEdgeRatios = mesh_->getCoedgeEdgeRatios();
  TEUCHOS_ASSERT( !coedgeEdgeRatios.is_null() );

  // loop over the local cells
  for ( int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation rel = (*cells[k]).relations();

      TEUCHOS_ASSERT_EQUALITY( 3, rel.size() );

      // extract the nodal coordinates
      Teuchos::Tuple<Point,3> localNodes = mesh_->getNodeCoordinates( rel );

      // loop over the edges
      for ( int l=0; l<3; l++ )
      {
          // global indices of subsequent nodes
          Teuchos::Tuple<int,2> nodeIndices;
          nodeIndices[0] = (*rel[ l       ].entity()).identifier() - 1;
          nodeIndices[1] = (*rel[ (l+1)%3 ].entity()).identifier() - 1;

          // co-edge / edge ratio
          double alpha = coedgeEdgeRatios[k][l];

          // -------------------------------------------------------------------
          // Compute the integral
          //
          //    I = \int_{x0}^{xj} (xj-x0).A(x) / |xj-x0| dx
          //
          // numerically by the midpoint rule, i.e.,
          //
          //    I ~ |xj-x0| * (xj-x0) . A( 0.5*(xj+x0) ) / |xj-x0|.
          //
          Point midpoint; // get A(midpoint)
          const Point & node0 = localNodes[ l       ];
          const Point & node1 = localNodes[ (l+1)%3 ];
          for (int i=0; i<midpoint.size(); i++ )
              midpoint[i] = 0.5 * ( node0[i] + node1[i] );

          // -------------------------------------------------------------------
          // Project vector field onto the edge.
          Teuchos::RCP<Point> a = mvp_->getA( midpoint );
          // Instead of first computing the projection over the normalized edge
          // and then multiply it with the edge length, don't normalize the
          // edge vector.
          double aInt = 0.0;
          for (int i=0; i<midpoint.size(); i++ )
              aInt += ( node1[i] - node0[i] ) * (*a)[i];

          // We'd like to insert the 2x2 matrix
          //
          //     [ alpha                     , - alpha * exp( -IM * aInt ) ]
          //     [ - alpha * exp( IM * aInt ), alpha                       ]
          //
          // at the indices   [ nodeIndices[0], nodeIndices[1] ] for every index pair
          // that shares and edge.
          // Do that now, just blockwise for real and imaginary part.
          Epetra_IntSerialDenseVector indices(4);
          indices[0] = 2*nodeIndices[0];
          indices[1] = 2*nodeIndices[0]+1;
          indices[2] = 2*nodeIndices[1];
          indices[3] = 2*nodeIndices[1]+1;

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
          // -------------------------------------------------------------------
      }
  }

  // calls FillComplete by default
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix.GlobalAssemble() );

  return;
}
// =============================================================================
const Epetra_FECrsGraph
Ginla::EpetraFVM::KeoFactory::
buildKeoGraph() const
{
  TEUCHOS_ASSERT( !mesh_.is_null() );
  Epetra_FECrsGraph keoGraph( Copy, *mesh_->getComplexMap(), 0 );

  std::vector<stk::mesh::Entity*> cells = mesh_->getOwnedCells();

  for ( int k=0; k<cells.size(); k++ )
  {
      // get the nodes local to the cell
      stk::mesh::PairIterRelation rel = (*cells[k]).relations();

      TEUCHOS_ASSERT_EQUALITY( 3, rel.size() );

      // extract the nodal coordinates
      Teuchos::Tuple<Point,3> localNodes = mesh_->getNodeCoordinates( rel );

      // loop over the edges
      for ( int l=0; l<3; l++ )
      {
          // global indices of subsequent nodes
          Teuchos::Tuple<int,2> nodeIndices;
          nodeIndices[0] = (*rel[ l       ].entity()).identifier() - 1;
          nodeIndices[1] = (*rel[ (l+1)%3 ].entity()).identifier() - 1;

          Epetra_IntSerialDenseVector indices(4);
          indices[0] = 2*nodeIndices[0];
          indices[1] = 2*nodeIndices[0]+1;
          indices[2] = 2*nodeIndices[1];
          indices[3] = 2*nodeIndices[1]+1;

          // sum it all in!
          TEUCHOS_ASSERT_EQUALITY( 0, keoGraph.InsertGlobalIndices( indices.Length(), indices.Values(),
                                                                    indices.Length(), indices.Values()
                                                                  )
                                 );
          // -------------------------------------------------------------------
      }
  }

  TEUCHOS_ASSERT_EQUALITY( 0, keoGraph.GlobalAssemble() );

  return keoGraph;
}
// =============================================================================
