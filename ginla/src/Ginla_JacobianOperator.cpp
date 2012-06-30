// @HEADER
//
//    Ginla Jacobian operator.
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
#include "Ginla_JacobianOperator.hpp"
#include "Ginla_KeoContainer.hpp"
#include "Ginla_StkMesh.hpp"
#include "Ginla_ScalarPotential_Virtual.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

namespace Ginla {
// =============================================================================
JacobianOperator::
JacobianOperator(const Teuchos::RCP<const Ginla::StkMesh> &mesh,
                 const Teuchos::RCP<Ginla::ScalarPotential::Virtual> &scalarPotential,
                 const double g,
                 const Teuchos::RCP<const Epetra_Vector> &thickness,
                 const Teuchos::RCP<Ginla::KeoContainer> &keoContainer,
                 const Teuchos::RCP<const Epetra_Vector> &current_X
                 ) :
  useTranspose_( false ),
  mesh_( mesh ),
  scalarPotential_( scalarPotential ),
  g_( g ),
  thickness_( thickness ),
  keoContainer_( keoContainer ),
  keo_( Teuchos::null ),
  current_X_( current_X ),
  diag0_( Teuchos::rcp( new Epetra_Vector( *(mesh->getComplexNonOverlapMap())) ) ),
  diag1b_( Teuchos::rcp( new Epetra_Vector( mesh->getControlVolumes()->Map()) ) )
{
  this->rebuildDiags_(*scalarPotential->get_p_init(), current_X_);
}
// =============================================================================
JacobianOperator::
~JacobianOperator ()
{
}
// =============================================================================
int
JacobianOperator::
SetUseTranspose( bool useTranspose )
{
  useTranspose_ = useTranspose;
  return 0;
}
// =============================================================================
int
JacobianOperator::
Apply(const Epetra_MultiVector &X,
      Epetra_MultiVector &Y
      ) const
{
  // Add the terms corresponding to the nonlinear terms.
  // A = K + I * thickness * (V + g * 2*|psi|^2)
  // B = g * diag( thickness * psi^2 )

#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keo_.is_null() );
#endif

  // K*psi
  TEUCHOS_ASSERT_EQUALITY(0, keo_->Apply(X, Y));

  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());
  int numMyPoints = controlVolumes.MyLength();

#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( 2*numMyPoints, X.MyLength() );
#endif

  for ( int vec=0; vec<X.NumVectors(); vec++ )
  {
    // For the parts Re(psi)Im(phi), Im(psi)Re(phi), the (2*k+1)th
    // component of X needs to be summed into the (2k)th component of Y,
    // likewise for (2k) -> (2k+1).
    // The Epetra class cannot currently handle this situation
    // (e.g., by Epetra_Vector::Multiply()), so we
    // need to access the vector entries one-by-one. And then, while
    // we're at it, let's include all the other terms in the loop
    // too. (It would actually be possible to have the terms
    // 2k/2k and 2k+1/2k+1 handled by Multiply().
    for ( int k=0; k<numMyPoints; k++ )
    {
      (*Y(vec))[2*k]   += (*diag0_) [2*k]   * X[vec][2*k]
                        + (*diag1b_)[k]     * X[vec][2*k+1];
      (*Y(vec))[2*k+1] += (*diag1b_)[k]     * X[vec][2*k]
                        + (*diag0_) [2*k+1] * X[vec][2*k+1];
    }

//        // add terms corresponding to  diag( psi^2 ) * \conj{phi}
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
  }

//    // take care of the shifting
//    if ( alpha_ != 0.0 || beta_ != -1.0 )
//    TEUCHOS_ASSERT_EQUALITY( 0, Y.Update( alpha_, X, -beta_ ) );

  return 0;
}
// =============================================================================
int
JacobianOperator::
ApplyInverse( const Epetra_MultiVector &X,
              Epetra_MultiVector &Y
              ) const
{
  TEUCHOS_TEST_FOR_EXCEPT( "Not implemented." );
  return -1;
}
// =============================================================================
double
JacobianOperator::
NormInf() const
{
  TEUCHOS_TEST_FOR_EXCEPT( "Not yet implemented." );
  return 0.0;
}
// =============================================================================
const char *
JacobianOperator::
Label() const
{
  return "Jacobian operator for Ginzburg--Landau";
}
// =============================================================================
bool
JacobianOperator::
UseTranspose() const
{
  return useTranspose_;
}
// =============================================================================
bool
JacobianOperator::
HasNormInf() const
{
  return false;
}
// =============================================================================
const Epetra_Comm &
JacobianOperator::
Comm() const
{
  return current_X_->Comm();
}
// =============================================================================
const Epetra_Map &
JacobianOperator::
OperatorDomainMap() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keo_.is_null() );
#endif
  return keo_->OperatorDomainMap();
}
// =============================================================================
const Epetra_Map &
JacobianOperator::
OperatorRangeMap() const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keo_.is_null() );
#endif
  return keo_->OperatorRangeMap();
}
// =============================================================================
void
JacobianOperator::
rebuild(const Teuchos::Array<double> &mvpParams,
        const Teuchos::Array<double> &spParams,
        const Teuchos::RCP<const Epetra_Vector> &current_X
        )
{
  // Rebuild diagonals immediately.
  this->rebuildDiags_(spParams, current_X);

  // Rebuild the KEO.
  keo_ = keoContainer_->getKeo(mvpParams);

  return;
}
// =============================================================================
void
JacobianOperator::
rebuildDiags_(const Teuchos::Array<double> &spParams,
              const Teuchos::RCP<const Epetra_Vector> &current_X
              )
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !current_X.is_null() );
  TEUCHOS_ASSERT( !scalarPotential_.is_null() );
#endif

  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());
  int numMyPoints = controlVolumes.MyLength();

  int kk;
  double val;
  int indices[2];
  double vals[2];
  for ( int k=0; k<numMyPoints; k++ )
  {
    // rebuild diag0
    double alpha = controlVolumes[k] * (*thickness_)[k] * (
        scalarPotential_->getV(k, spParams)
        + g_ * 2.0 *
          ( (*current_X)[2*k] * (*current_X)[2*k]
          + (*current_X)[2*k+1]*(*current_X)[2*k+1] )
        );
    double realX2 = g_ * controlVolumes[k] * (*thickness_)[k] * (
          (*current_X)[2*k]  *(*current_X)[2*k]
        - (*current_X)[2*k+1]*(*current_X)[2*k+1]
        );
    vals[0]    = alpha + realX2;
    vals[1]    = alpha - realX2;
    indices[0] = 2*k;
    indices[1] = 2*k+1;
    TEUCHOS_ASSERT_EQUALITY( 0, diag0_->ReplaceMyValues( 2, vals, indices ) );

    // rebuild diag1b
    double imagX2 = g_ * controlVolumes[k] * (*thickness_)[k] * (
        2.0 * (*current_X)[2*k] * (*current_X)[2*k+1]
        );
    TEUCHOS_ASSERT_EQUALITY( 0, diag1b_->ReplaceMyValues( 1, &imagX2, &k ) );
  }

  current_X_ = current_X;
  return;
}
// =============================================================================
} // namespace Ginla
