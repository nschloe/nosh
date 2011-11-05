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
#include "Ginla_JacobianOperator.hpp"
#include "Ginla_KeoFactory.hpp"
#include "Ginla_StkMesh.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Epetra_Vector.h>

namespace Ginla {
// =============================================================================
JacobianOperator::
JacobianOperator( const Teuchos::RCP<Ginla::StkMesh>                 & mesh,
                  const Teuchos::RCP<const Epetra_Vector>            & thickness,
                  const Teuchos::RCP<Ginla::MagneticVectorPotential> & mvp,
                  const Teuchos::RCP<Epetra_Vector>                  & current_X
                ):
        useTranspose_( false ),
        comm_( mesh->getComm() ),
        mesh_( mesh ),
        thickness_( thickness ),
        keoFactory_( Teuchos::rcp( new Ginla::KeoFactory( mesh, thickness, mvp ) ) ),
        keoMatrix_( Teuchos::rcp( new Epetra_FECrsMatrix( Copy, keoFactory_->buildKeoGraph() ) ) ),
        current_X_ ( current_X ),
        temperature_( 0.0 ),
        isDiag0UpToDate_( false ),
        diag0_ ( Teuchos::rcp( new Epetra_Vector(*mesh->getComplexNonOverlapMap()) ) ),
        isDiag1UpToDate_( false ),
        diag1a_( Teuchos::rcp( new Epetra_Vector(mesh->getControlVolumes()->Map()) ) ),
        diag1b_( Teuchos::rcp( new Epetra_Vector(mesh->getControlVolumes()->Map()) ) )
{
    // Fill the matrix immediately.
    keoFactory_->buildKeo( *keoMatrix_ );
}
// =============================================================================
JacobianOperator::
~JacobianOperator ()
{
}
// =============================================================================
int
JacobianOperator::
SetUseTranspose( bool UseTranspose )
{
    useTranspose_ = UseTranspose;
    return 0;
}
// =============================================================================
int
JacobianOperator::
Apply ( const Epetra_MultiVector & X,
              Epetra_MultiVector & Y
      ) const
{
    // Add the terms corresponding to the nonlinear terms.
    // A = K - I * thickness * ( (1-temp) - 2*|psi|^2 )
    // B = diag( thickness * psi^2 )

    TEUCHOS_ASSERT( !keoMatrix_.is_null() );
    TEUCHOS_ASSERT( !current_X_.is_null() );

    // K*psi
    TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix_->Apply( X, Y ) );

    const Epetra_Vector & controlVolumes =  *(mesh_->getControlVolumes());
    int numMyPoints = controlVolumes.MyLength();

    TEUCHOS_ASSERT_EQUALITY( 2*numMyPoints, X.MyLength() );

    // add terms corresponding to (1-T-2*|psi|^2) * phi
    if ( !isDiag0UpToDate_ )
        this->rebuildDiag0_();
    for ( int vec=0; vec<X.NumVectors(); vec++ )
        TEUCHOS_ASSERT_EQUALITY( 0, Y(vec)->Multiply( 1.0, *diag0_, *(X(vec)), 1.0 ) );

    // cache rePhiSquare, imPhiSquare
    if ( !isDiag1UpToDate_ )
        this->rebuildDiag1_();

    // Add terms corresponding to  diag( psi^2 ) * \conj{phi}.
    // There's no (visible) speedup when compared to computing the coefficients in every step,
    // but there is the cost of two additional vectors. Think about going back to the old
    // setup (below).
    for ( int vec=0; vec<X.NumVectors(); vec++ )
    {
        for ( int k=0; k<numMyPoints; k++ )
        {
            // real part
            TEUCHOS_ASSERT_EQUALITY( 0, Y.SumIntoMyValue( 2*k,   vec, - (*diag1a_)[k] * X[vec][2*k] - (*diag1b_)[k] * X[vec][2*k+1] ) );
            // imaginary part
            TEUCHOS_ASSERT_EQUALITY( 0, Y.SumIntoMyValue( 2*k+1, vec, - (*diag1b_)[k] * X[vec][2*k] + (*diag1a_)[k] * X[vec][2*k+1] ) );
            // There is virtually no difference in speed between calling one or
            // two SumIntoMyValue()s per iteration. Hence, prefer the version
            // with less code. The one-call version is kept for reference.
//            int indices[2];
//            double values[2];
//            indices[0] = 2*k;
//            indices[1] = 2*k+1;
//            values[0] = - (*diag1a_)[k] * X[vec][2*k] - (*diag1b_)[k] * X[vec][2*k+1];
//            values[1] = - (*diag1b_)[k] * X[vec][2*k] + (*diag1a_)[k] * X[vec][2*k+1];
//            TEUCHOS_ASSERT_EQUALITY( 0, Y(vec)->SumIntoMyValues( 2, values, indices )  );
        }
    }

//    // add terms corresponding to  diag( psi^2 ) * \conj{phi}
//    for ( int vec=0; vec<X.NumVectors(); vec++ )
//    {
//        for ( int k=0; k<numMyPoints; k++ )
//        {
//            double rePhiSquare = controlVolumes[k] * (*thickness_)[k] * (
//                                 (*current_X_)[2*k]*(*current_X_)[2*k] - (*current_X_)[2*k+1]*(*current_X_)[2*k+1]
//                                 );
//            // Im(phi^2)
//            double imPhiSquare = controlVolumes[k] * (*thickness_)[k] * (
//                                 2.0 * (*current_X_)[2*k] * (*current_X_)[2*k+1]
//                                );
//            // real part
//            TEUCHOS_ASSERT_EQUALITY( 0, Y.SumIntoMyValue( 2*k,   vec, - rePhiSquare * X[vec][2*k] - imPhiSquare * X[vec][2*k+1] ) );
//            // imaginary part
//            TEUCHOS_ASSERT_EQUALITY( 0, Y.SumIntoMyValue( 2*k+1, vec, - imPhiSquare * X[vec][2*k] + rePhiSquare * X[vec][2*k+1] ) );
//        }
//    }

//    // take care of the shifting
//    if ( alpha_ != 0.0 || beta_ != -1.0 )
//    TEUCHOS_ASSERT_EQUALITY( 0, Y.Update( alpha_, X, -beta_ ) );

    return 0;
}
// =============================================================================
int
JacobianOperator::
ApplyInverse ( const Epetra_MultiVector & X,
                     Epetra_MultiVector & Y
             ) const
{
    TEUCHOS_TEST_FOR_EXCEPT( "Not implemented." );
    return -1;
}
// =============================================================================
double
JacobianOperator::
NormInf () const
{
    TEUCHOS_TEST_FOR_EXCEPT( "Not yet implemented." );
    return 0.0;
}
// =============================================================================
const char *
JacobianOperator::
Label () const
{
    return "Jacobian operator for Ginzburg--Landau";
}
// =============================================================================
bool
JacobianOperator::
UseTranspose () const
{
    return useTranspose_;
}
// =============================================================================
bool
JacobianOperator::
HasNormInf () const
{
    return false;
}
// =============================================================================
const Epetra_Comm &
JacobianOperator::
Comm () const
{
    return comm_;
}
// =============================================================================
const Epetra_Map &
JacobianOperator::
OperatorDomainMap () const
{
    TEUCHOS_ASSERT( !keoMatrix_.is_null() );
    return keoMatrix_->OperatorDomainMap();
}
// =============================================================================
const Epetra_Map &
JacobianOperator::
OperatorRangeMap () const
{
    TEUCHOS_ASSERT( !keoMatrix_.is_null() );
    return keoMatrix_->OperatorRangeMap();
}
// =============================================================================
void
JacobianOperator::
rebuild( const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams,
         const double temperature,
         const Teuchos::RCP<const Epetra_Vector> & current_X
       )
{
    // set the entities for the operator
    temperature_ = temperature;
    current_X_ = current_X;
    isDiag0UpToDate_ = false;
    isDiag1UpToDate_ = false;

    // rebuild the keo
    keoFactory_->updateParameters( mvpParams );
    keoFactory_->buildKeo( *keoMatrix_ );

    return;
}
// =============================================================================
void
JacobianOperator::
rebuildDiag0_() const
{
    const Epetra_Vector & controlVolumes = *(mesh_->getControlVolumes());
    int numMyPoints = controlVolumes.MyLength();

    // build scaling vector
    for ( int k=0; k<numMyPoints; k++ )
    {
        int numEntries = 2;
        double alpha = controlVolumes[k] * (*thickness_)[k] * (
                       1.0 - temperature_
                       - 2.0 * ( (*current_X_)[2*k]*(*current_X_)[2*k] + (*current_X_)[2*k+1]*(*current_X_)[2*k+1] )
                       );
        double vals[2];
        vals[0] = alpha; vals[1] = alpha;
        int indices[2];
        indices[0] = 2*k; indices[1] = 2*k+1;
        TEUCHOS_ASSERT_EQUALITY( 0, diag0_->ReplaceMyValues( numEntries, vals, indices ) );
    }
    isDiag0UpToDate_ = true;

    return;
}
// =============================================================================
void
JacobianOperator::
rebuildDiag1_() const
{
    const Epetra_Vector & controlVolumes = *(mesh_->getControlVolumes());
    int numMyPoints = controlVolumes.MyLength();

    for ( int k=0; k<numMyPoints; k++ )
    {
        double rePhiSquare = controlVolumes[k] * (*thickness_)[k] * (
                             (*current_X_)[2*k]*(*current_X_)[2*k] - (*current_X_)[2*k+1]*(*current_X_)[2*k+1]
                             );
        TEUCHOS_ASSERT_EQUALITY( 0, diag1a_->ReplaceMyValues( 1, &rePhiSquare, &k ) );
        // Im(phi^2)
        double imPhiSquare = controlVolumes[k] * (*thickness_)[k] * (
                             2.0 * (*current_X_)[2*k] * (*current_X_)[2*k+1]
                             );
        TEUCHOS_ASSERT_EQUALITY( 0, diag1b_->ReplaceMyValues( 1, &imPhiSquare, &k ) );
    }
    isDiag1UpToDate_ = true;
    return;
}
// =============================================================================
} // namespace Ginla
