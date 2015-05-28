// @HEADER
//
//    Nosh Jacobian operator.
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
#include "nosh/JacobianOperator.hpp"

#include <map>
#include <string>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "nosh/ParameterMatrix_Virtual.hpp"
#include "nosh/StkMesh.hpp"
#include "nosh/ScalarField_Virtual.hpp"

namespace Nosh
{
// =============================================================================
JacobianOperator::
JacobianOperator(
    const std::shared_ptr<const Nosh::StkMesh> &mesh,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &scalarPotential,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
    const std::shared_ptr<Nosh::ParameterMatrix::Virtual> &matrix
    ) :
  useTranspose_(false),
  mesh_(mesh),
  scalarPotential_(scalarPotential),
  thickness_(thickness),
  keo_(matrix),
  diag0_(*(mesh->getComplexNonOverlapMap())),
  diag1b_(mesh->getControlVolumes()->Map())
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
SetUseTranspose(bool useTranspose)
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
  // B = g * diag(thickness * psi^2)

  const int numMyPoints = mesh_->getControlVolumes()->MyLength();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, X.MyLength());
#endif

  // Y = K*X
  TEUCHOS_ASSERT_EQUALITY(0, keo_->Apply(X, Y));

  for (int vec = 0; vec < X.NumVectors(); vec++) {
    // For the parts Re(psi)Im(phi), Im(psi)Re(phi), the (2*k+1)th
    // component of X needs to be summed into the (2k)th component of Y,
    // likewise for (2k) -> (2k+1).
    // The Epetra class cannot currently handle this situation
    // (e.g., by Epetra_Vector::Multiply()), so we
    // need to access the vector entries one-by-one. And then, while
    // we're at it, let's include all the other terms in the loop
    // too. (It would actually be possible to have the terms
    // 2k/2k and 2k+1/2k+1 handled by Multiply().
    for (int k = 0; k < numMyPoints; k++) {
      (*Y(vec))[2*k]  += diag0_ [2*k]   * X[vec][2*k]
                       + diag1b_[k]     * X[vec][2*k+1];
      (*Y(vec))[2*k+1] += diag1b_[k]     * X[vec][2*k]
                        + diag0_ [2*k+1] * X[vec][2*k+1];
    }
  }

//    // take care of the shifting
//    if (alpha_ != 0.0 || beta_ != -1.0)
//    TEUCHOS_ASSERT_EQUALITY(0, Y.Update(alpha_, X, -beta_));

  return 0;
}
// =============================================================================
int
JacobianOperator::
ApplyInverse(const Epetra_MultiVector &X,
              Epetra_MultiVector &Y
           ) const
{
  (void) X;
  (void) Y;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
}
// =============================================================================
double
JacobianOperator::
NormInf() const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not yet implemented.");
}
// =============================================================================
const char *
JacobianOperator::
Label() const
{
  return "Jacobian operator for nonlinear Schr\"odinger";
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
  return keo_->Comm();
}
// =============================================================================
const Epetra_Map &
JacobianOperator::
OperatorDomainMap() const
{
  return keo_->OperatorDomainMap();
}
// =============================================================================
const Epetra_Map &
JacobianOperator::
OperatorRangeMap() const
{
  return keo_->OperatorRangeMap();
}
// =============================================================================
void
JacobianOperator::
rebuild(
    const std::map<std::string, double> params,
    const Epetra_Vector & current_X
    )
{
  // Fill the KEO.
  // It is certainly a debatable design decision to have our own KEO in
  // JacobianOperator and not live of the cache of the builder. On the one
  // hand, in a typical continuation context, the same matrix is used in
  // computeF, the preconditioner, and here. It would then be sufficient to
  // store the matrix at one common place (e.g., the builder).  This might
  // however lead to complications in a situation like the following:
  //
  //   1. The matrix is requested by the Jacobian operator,
  //      a pointer to the common storage place is
  //      returned.
  //   2. The matrix is requested by computeF()
  //      with different parameters.
  //      Now also the Jacobian operato's instance
  //      has the altered matrix.
  //   3. The Jacobian operator uses the matrix.
  //
  // One might argue that this situation is does not occur in the given
  // context, but really the code shouldn't make any assumptions about it.
  // Besides, the matrix copy that happens in fill is not of much concern
  // computationally. Should this ever become an issue, revisit.
  keo_->setParameters(params);

  // Rebuild diagonals.
  this->rebuildDiags_(params, current_X);

  return;
}
// =============================================================================
void
JacobianOperator::
rebuildDiags_(
    const std::map<std::string, double> params,
    const Epetra_Vector &x
    )
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(scalarPotential_);
#endif

  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());

  const double g = params.at("g");

  const Epetra_Vector thicknessValues = thickness_->getV(params);
#ifndef NDEBUG
  TEUCHOS_ASSERT(controlVolumes.Map().SameAs(thicknessValues.Map()));
#endif

  const Epetra_Vector scalarPotentialValues = scalarPotential_->getV(params);
#ifndef NDEBUG
  TEUCHOS_ASSERT(controlVolumes.Map().SameAs(scalarPotentialValues.Map()));
#endif

  for (int k = 0; k < controlVolumes.MyLength(); k++) {
    // rebuild diag0
    const double alpha = controlVolumes[k] * thicknessValues[k]
                         * (scalarPotentialValues[k]
                             + g * 2.0 * (x[2*k]*x[2*k] + x[2*k+1]*x[2*k+1])
                          );
    const double realX2 = g * controlVolumes[k] * thicknessValues[k]
                          * (x[2*k]*x[2*k] - x[2*k+1]*x[2*k+1]);
    diag0_[2*k]   = alpha + realX2;
    diag0_[2*k+1] = alpha - realX2;

    // rebuild diag1b
    diag1b_[k] = g * controlVolumes[k] * thicknessValues[k]
                 * (2.0 * x[2*k] * x[2*k+1]);
  }

  return;
}
// =============================================================================
} // namespace Nosh
