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

#include "Ginla_FVM_KineticEnergyOperator.h"

// =============================================================================
Ginla::FVM::KineticEnergyOperator::
KineticEnergyOperator( const Teuchos::RCP<VIO::Mesh::Mesh>                         & mesh,
                       const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp
                     ):
        mesh_ ( mesh ),
        mvp_( mvp ),
        keoGraph_( Teuchos::null ),
        keo_( Teuchos::null ),
        mu_ ( 0.0 ),
        scaling_( Teuchos::tuple( 0.0, 0.0, 0.0 ) ),
        keoMu_( 0.0 ),
        keoScaling_( Teuchos::tuple( 0.0, 0.0, 0.0 ) )
{
}
// =============================================================================
Ginla::FVM::KineticEnergyOperator::
~KineticEnergyOperator()
{
}
// =============================================================================
int
Ginla::FVM::KineticEnergyOperator::
SetUseTranspose( bool UseTranspose )
{
    useTranspose_ = UseTranspose;
    return 0;
}
// =============================================================================
int
Ginla::FVM::KineticEnergyOperator::
Apply ( const Epetra_MultiVector & X,
              Epetra_MultiVector & Y
      ) const
{
    // reassemble linear operator if necessary
    if ( !this->keoUpToDate_() )
        this->assembleKeo_();

    TEUCHOS_ASSERT_EQUALITY( 0, keo_->Apply( X, Y ) );

    return 0;
}
// =============================================================================
int
Ginla::FVM::KineticEnergyOperator::
ApplyInverse ( const Epetra_MultiVector & X,
                     Epetra_MultiVector & Y
             ) const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic,
                        "Not yet implemented." );
    return -1;
}
// =============================================================================
double
Ginla::FVM::KineticEnergyOperator::
NormInf () const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic,
                        "Not yet implemented." );
    return 0.0;
}
// =============================================================================
const char *
Ginla::FVM::KineticEnergyOperator::
Label () const
{
    return "Kinetic energy operator for Ginzburg--Landau";
}
// =============================================================================
bool
Ginla::FVM::KineticEnergyOperator::
UseTranspose () const
{
    return useTranspose_;
}
// =============================================================================
bool
Ginla::FVM::KineticEnergyOperator::
HasNormInf () const
{
    return false;
}
// =============================================================================
const Epetra_Comm &
Ginla::FVM::KineticEnergyOperator::
Comm () const
{
    return comm_;
}
// =============================================================================
const Epetra_Map &
Ginla::FVM::KineticEnergyOperator::
OperatorDomainMap () const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic,
                        "Not yet implemented." );
}
// =============================================================================
const Epetra_Map &
Ginla::FVM::KineticEnergyOperator::
OperatorRangeMap () const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic,
                        "Not yet implemented." );
}
// =============================================================================
void
Ginla::FVM::KineticEnergyOperator::
setParameters( const double mu,
               const Teuchos::Tuple<double,3> & scaling
             )
{
    mu_ = mu;
    scaling_ = scaling;
    return;
}
// =============================================================================
void
Ginla::FVM::KineticEnergyOperator::
assembleKeo_() const
{
  if ( keoGraph_.is_null() )
      this->createKeoGraph_();

  keo_  = Teuchos::rcp( new Epetra_FECrsMatrix( Copy, *keoGraph_ ) );
  keo_->PutScalar( 0.0 );

  // set scaling and external magnetic field
  mesh_->scale( scaling_ );
  mvp_->setMu( mu_ );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Loop over the elements, create local load vector and mass matrix,
  // and insert them into the global matrix.
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems = mesh_->getElems();
  Teuchos::ArrayRCP<Point> nodes = mesh_->getNodes();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeLengths = mesh_->getEdgeLengths();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > coedgeLengths = mesh_->getCoedgeLengths();

  TEUCHOS_ASSERT( !elems.is_null() );
  TEUCHOS_ASSERT( !nodes.is_null() );
  TEUCHOS_ASSERT( !edgeLengths.is_null() );
  TEUCHOS_ASSERT( !coedgeLengths.is_null() );

  // loop over the local elements
  for ( int k=0; k<elems.size(); k++ )
  {
      Teuchos::ArrayRCP<ORD> & elem = elems[k];

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
          //    I = \int_{x0}^{xj} (xj-x0).A(x) dx
          //
          // numerically by the midpoint rule, i.e.,
          //
          //    I ~ |xj-x0| * (xj-x0) . A( 0.5*(xj+x0) ).
          //
          Point midpoint; // get A(midpoint)
          Point & node0 = nodes[ nodeIndices[0] ];
          Point & node1 = nodes[ nodeIndices[1] ];
          for (int i=0; i<midpoint.size(); i++ )
              midpoint[i] = 0.5 * ( node0[i] + node1[i] );

          // -------------------------------------------------------------------
          // Project vector field onto the edge.
          Teuchos::RCP<Point> a = mvp_->getA( midpoint );
          double aProj = 0.0;
          for (int i=0; i<midpoint.size(); i++ )
              aProj += ( node1[i] - node0[i] ) * (*a)[i];

          const double aInt = aProj * edgeLengths[k][l];

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
          values(0,0) = alpha;
          values(0,1) = 0.0;
          values(1,0) = 0.0;
          values(1,1) = alpha;

          double alphaCosAInt = alpha * cos(aInt);
          double alphaSinAInt = alpha * sin(aInt);

          values(0,2) = - alphaCosAInt;
          values(0,3) = - alphaSinAInt;
          values(1,2) =   alphaSinAInt;
          values(1,3) = - alphaCosAInt;

          values(2,0) = - alphaCosAInt;
          values(2,1) =   alphaSinAInt;
          values(3,0) = - alphaSinAInt;
          values(3,1) = - alphaCosAInt;

          values(2,2) = alpha;
          values(2,3) = 0.0;
          values(3,2) = 0.0;
          values(3,3) = alpha;

          // sum it all in!
          keo_->SumIntoGlobalValues ( indices, values );
          // -------------------------------------------------------------------
      }
  }

  // TODO both necessary?
  keo_->GlobalAssemble();
  keo_->FillComplete();

  keoMu_ = mu;
  keoScaling_ = scaling;

  return;
}
// =============================================================================
void
Ginla::FVM::KineticEnergyOperator::
createKeoGraph_() const
{
  keoGraph_ = Teuchos::rcp( new Epetra_FECrsGraph( Copy, *(komplex_->getRealMap()), 0 ) );

  TEUCHOS_ASSERT( !mesh_.is_null() );
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems = mesh_->getElems();

  for ( int k=0; k<elems.size(); k++ )
  {
      Teuchos::ArrayRCP<ORD> & elem = elems[k];

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
          keoGraph_->InsertGlobalIndices( indices.Length(), indices.Values(),
                                          indices.Length(), indices.Values()
                                        );
          // -------------------------------------------------------------------
      }
  }

  keoGraph_->FillComplete();

  return;
}
// =============================================================================
bool
Ginla::FVM::KineticEnergyOperator::
keoUpToDate_() const
{
    return    !keo_.is_null()
           && keoMu_         == mu_
           && keoScaling_[0] == scaling_[0]
           && keoScaling_[1] == scaling_[1]
           && keoScaling_[2] == scaling_[2];
}
// =============================================================================
