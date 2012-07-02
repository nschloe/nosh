// @HEADER
//
//    Cuantico Jacobian operator.
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
#include "Cuantico_JacobianOperator.hpp"
#include "Cuantico_KeoContainer.hpp"
#include "Cuantico_StkMesh.hpp"
#include "Cuantico_ScalarPotential_Virtual.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

namespace Cuantico {
// =============================================================================
JacobianOperator::
JacobianOperator(const Teuchos::RCP<const Cuantico::StkMesh> &mesh,
                 const Teuchos::RCP<const Cuantico::ScalarPotential::Virtual> &scalarPotential,
                 const Teuchos::RCP<const Epetra_Vector> &thickness,
                 const Teuchos::RCP<const Cuantico::KeoContainer> &keoContainer
                 ) :
  useTranspose_( false ),
  mesh_( mesh ),
  scalarPotential_( scalarPotential ),
  thickness_( thickness ),
  keoContainer_( keoContainer ),
  keo_( Teuchos::null ),
  diag0_(Teuchos::rcp(new Epetra_Vector(*(mesh->getComplexNonOverlapMap())))),
  diag1b_(Teuchos::rcp(new Epetra_Vector(mesh->getControlVolumes()->Map())))
{
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
  return "Jacobian operator for nonlinear Schrödinger";
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
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !keo_.is_null() );
#endif
  return keo_->Comm();
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
rebuild(const double g,
        const Teuchos::Array<double> &spParams,
        const Teuchos::Array<double> &mvpParams,
        const Teuchos::RCP<const Epetra_Vector> &current_X
        )
{
  // Rebuild the KEO.
  keo_ = keoContainer_->getKeo(mvpParams);

  // Rebuild diagonals.
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !current_X.is_null() );
#endif
  this->rebuildDiags_(g, spParams, *current_X);

  return;
}
// =============================================================================
void
JacobianOperator::
rebuildDiags_(const double g,
              const Teuchos::Array<double> &spParams,
              const Epetra_Vector &x
              )
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !scalarPotential_.is_null() );
#endif

  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());

  for (int k=0; k<controlVolumes.MyLength(); k++)
  {
    // rebuild diag0
    const double alpha = controlVolumes[k] * (*thickness_)[k]
                       * ( scalarPotential_->getV(k, spParams)
                         + g * 2.0 * (x[2*k]*x[2*k] + x[2*k+1]*x[2*k+1])
                         );
    const double realX2 = g * controlVolumes[k] * (*thickness_)[k]
                        * ( x[2*k]*x[2*k] - x[2*k+1]*x[2*k+1] );
    (*diag0_)[2*k]   = alpha + realX2;
    (*diag0_)[2*k+1] = alpha - realX2;

    // rebuild diag1b
    (*diag1b_)[k] = g * controlVolumes[k] * (*thickness_)[k]
                  * (2.0 * x[2*k] * x[2*k+1]);
  }

  return;
}
// =============================================================================
} // namespace Cuantico
