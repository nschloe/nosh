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

#include "Ginla_EpetraFVM_KeoFactory.h"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Comm.h>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosBlockCGSolMgr.hpp>

#include <ml_epetra_preconditioner.h>

// =============================================================================
// some typdefs for Belos
typedef double                           ST;
typedef Epetra_MultiVector               MV;
typedef Epetra_Operator                  OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
// =============================================================================
Ginla::EpetraFVM::KeoFactory::
KeoFactory( const Teuchos::RCP<VIO::EpetraMesh::Mesh>                   & mesh,
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
buildKeo( Epetra_FECrsMatrix & keoMatrix,
          const double mu,
          const Teuchos::Tuple<double,3> & scaling
        ) const
{
  keoMatrix.PutScalar( 0.0 );

  // set teh parameters
  mvp_->setMu( mu );
  mesh_->scale( scaling );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the elements, create local load vector and mass matrix,
  // and insert them into the global matrix.
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> > elems = mesh_->getElems();
  Teuchos::ArrayRCP<const Point> nodes = mesh_->getNodes();

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeLengths = mesh_->getEdgeLengths();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > coedgeLengths = mesh_->getCoedgeLengths();

  TEUCHOS_ASSERT( !elems.is_null() );
  TEUCHOS_ASSERT( !nodes.is_null() );
  TEUCHOS_ASSERT( !edgeLengths.is_null() );
  TEUCHOS_ASSERT( !coedgeLengths.is_null() );

  // loop over the local elements
  for ( int k=0; k<elems.size(); k++ )
  {
      Teuchos::ArrayRCP<int> & elem = elems[k];

      // loop over the edges
      int n = elem.size();
      for ( int l=0; l<n; l++ )
      {

          // two subsequent indices determine an edge
          Teuchos::Tuple<ORD,2> nodeIndices;
          nodeIndices[0] = elem[l];
          nodeIndices[1] = elem[(l+1)%n];

          // co-edge / edge ratio
          double alpha = coedgeLengths[k][l]
                       / edgeLengths[k][l];

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
          const Point & node0 = nodes[ nodeIndices[0] ];
          const Point & node1 = nodes[ nodeIndices[1] ];
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
  Epetra_FECrsGraph keoGraph( Copy, *(mesh_->getComplexValuesMap()), 0 );
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> > elems = mesh_->getElems();

  for ( int k=0; k<elems.size(); k++ )
  {
      Teuchos::ArrayRCP<int> & elem = elems[k];

      // loop over the edges
      int n = elem.size();
      for ( int l=0; l<n; l++ )
      {
          // two subsequent indices determine an edge
          Teuchos::Tuple<ORD,2> nodeIndices;
          nodeIndices[0] = elem[l];
          nodeIndices[1] = elem[(l+1)%n];

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
