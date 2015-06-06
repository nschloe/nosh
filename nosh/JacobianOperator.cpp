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

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include "nosh/ParameterMatrix_Keo.hpp"
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
    const std::shared_ptr<Nosh::ParameterMatrix::Keo> &keo
    ) :
  mesh_(mesh),
  scalarPotential_(scalarPotential),
  thickness_(thickness),
  keo_(keo),
  diag0_(Teuchos::rcp(mesh->getComplexNonOverlapMap())),
  diag1b_(mesh->getControlVolumes()->getMap())
{
}
// =============================================================================
JacobianOperator::
~JacobianOperator ()
{
}
// =============================================================================
void
JacobianOperator::
apply(
    const Tpetra::MultiVector<double,int,int> &X,
    Tpetra::MultiVector<double,int,int> &Y,
    Teuchos::ETransp mode,
    double alpha,
    double beta
    ) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      mode != Teuchos::NO_TRANS,
      "Only untransposed applies supported."
      );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      alpha != 1.0,
      "Only alpha==1.0 supported."
      );
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      beta != 0.0,
      "Only beta==0.0 supported."
      )
  // Add the terms corresponding to the nonlinear terms.
  // A = K + I * thickness * (V + g * 2*|psi|^2)
  // B = g * diag(thickness * psi^2)

  // Y = K*X
  keo_->apply(X, Y);

  const std::size_t numMyPoints = mesh_->getControlVolumes()->getLocalLength();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, X.getLocalLength());
#endif

  auto d0Data = diag0_.getData();
  auto d1bData = diag1b_.getData();

#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, d0Data.size());
  TEUCHOS_ASSERT_EQUALITY(numMyPoints, d1bData.size());
#endif

  for (std::size_t i = 0; i < X.getNumVectors(); i++) {
    auto xData = X.getVector(i)->getData();
    auto yData = Y.getVectorNonConst(i)->getDataNonConst();
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, xData.size());
    TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, yData.size());
#endif
    // For the parts Re(psi)Im(phi), Im(psi)Re(phi), the (2*k+1)th component of
    // X needs to be summed into the (2k)th component of Y, likewise for (2k)
    // -> (2k+1).
    // The Epetra class cannot currently handle this situation (e.g., by
    // Tpetra::Vector<double,int,int>::Multiply()), so we need to access the
    // vector entries one-by-one. And then, while we're at it, let's include
    // all the other terms in the loop too. (It would actually be possible to
    // have the terms 2k/2k and 2k+1/2k+1 handled by Multiply().
    for (std::size_t k = 0; k < numMyPoints; k++) {
      yData[2*k] += d0Data[2*k] * xData[2*k]
                  + d1bData[k]  * xData[2*k+1];
      yData[2*k+1] += d1bData[k]    * xData[2*k]
                    + d0Data[2*k+1] * xData[2*k+1];
    }
  }

//    // take care of the shifting
//    if (alpha_ != 0.0 || beta_ != -1.0)
//    TEUCHOS_ASSERT_EQUALITY(0, Y.Update(alpha_, X, -beta_));

  return;
}
// =============================================================================
Teuchos::RCP<const Tpetra::Map<int,int>>
JacobianOperator::
getDomainMap() const
{
  return keo_->getDomainMap();
}
// =============================================================================
Teuchos::RCP<const Tpetra::Map<int,int>>
JacobianOperator::
getRangeMap() const
{
  return keo_->getRangeMap();
}
// =============================================================================
void
JacobianOperator::
rebuild(
    const std::map<std::string, double> params,
    const Tpetra::Vector<double,int,int> & current_X
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
    const Tpetra::Vector<double,int,int> &x
    )
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(scalarPotential_);
#endif

  const auto & controlVolumes = *(mesh_->getControlVolumes());

  const double g = params.at("g");

  const auto thicknessValues = thickness_->getV(params);
#ifndef NDEBUG
  TEUCHOS_ASSERT(controlVolumes.getMap()->isSameAs(*thicknessValues.getMap()));
#endif

  const auto scalarPotentialValues = scalarPotential_->getV(params);
#ifndef NDEBUG
  TEUCHOS_ASSERT(
      controlVolumes.getMap()->isSameAs(*scalarPotentialValues.getMap())
      );
#endif

  auto xData = x.getData();
  auto cData = controlVolumes.getData();
  auto tData = thicknessValues.getData();
  auto sData = scalarPotentialValues.getData();

  auto d0Data = diag0_.getDataNonConst();
  auto d1bData = diag1b_.getDataNonConst();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(cData.size(), tData.size());
  TEUCHOS_ASSERT_EQUALITY(tData.size(), sData.size());
  TEUCHOS_ASSERT_EQUALITY(2*sData.size(), xData.size());
  TEUCHOS_ASSERT_EQUALITY(sData.size(), d1bData.size());
#endif

  for (decltype(cData)::size_type k = 0; k < cData.size(); k++) {
    // rebuild diag0
    const double alpha = cData[k] * tData[k]
      * (sData[k]
          + g * 2.0 * (xData[2*k]*xData[2*k] + xData[2*k+1]*xData[2*k+1])
        );
    const double realX2 = g * cData[k] * tData[k]
      * (xData[2*k]*xData[2*k] - xData[2*k+1]*xData[2*k+1]);
    d0Data[2*k]   = alpha + realX2;
    d0Data[2*k+1] = alpha - realX2;

    // rebuild diag1b
    d1bData[k] = g * cData[k] * tData[k] * (2.0 * xData[2*k] * xData[2*k+1]);
  }

  return;
}
// =============================================================================
} // namespace Nosh
