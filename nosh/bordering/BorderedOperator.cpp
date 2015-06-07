// @HEADER
//
//    Regularized kinetic energy operator.
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

#include "BorderedOperator.hpp"

#include <vector>

#include "BorderingHelpers.hpp"

// =============================================================================
namespace Nosh
{
// =============================================================================
BorderedOperator::
BorderedOperator(const std::shared_ptr<Tpetra::Operator<double,int,int>> & innerOperator,
                 const Tpetra::Vector<double,int,int> & b,
                 const Tpetra::Vector<double,int,int> & c,
                 const double d
                ):
  innerOperator_(innerOperator),
  b_(b),
  c_(c),
  d_(d),
  useTranspose_(false),
  domainMap_(*Nosh::BorderingHelpers::extendMapBy1(innerOperator_->getDomainMap())),
  rangeMap_(*Nosh::BorderingHelpers::extendMapBy1(innerOperator_->getRangeMap()))
{
}
// =============================================================================
BorderedOperator::
~BorderedOperator()
{
}
// =============================================================================
int
BorderedOperator::
SetUseTranspose(bool useTranspose)
{
  innerOperator_->SetUseTranspose(useTranspose);
  useTranspose_ = useTranspose;
  return 0;
}
// =============================================================================
int
BorderedOperator::
Apply(const Tpetra::MultiVector<double,int,int> &X,
      Tpetra::MultiVector<double,int,int> &Y
     ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(X.getMap().SameAs(domainMap_));
  TEUCHOS_ASSERT(Y.getMap().SameAs(rangeMap_));
#endif
  const int n = X.NumVectors();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(n, Y.NumVectors());
#endif
  // Dissect X.
  Tpetra::Vector<double,int,int> innerX(innerOperator_->getDomainMap());
  std::vector<double> lambda(n);
  Nosh::BorderingHelpers::dissect(X, innerX, &lambda[0]);
  // Apply inner operator.
  Tpetra::Vector<double,int,int> innerY(innerOperator_->getRangeMap());
  TEUCHOS_ASSERT_EQUALITY(0, innerOperator_->Apply(innerX, innerY));

  const Tpetra::Vector<double,int,int> & rightBordering = useTranspose_? c_ : b_;
  const Tpetra::Vector<double,int,int> & lowerBordering = useTranspose_? b_ : c_;

  // Add right bordering.
  for (int k = 0; k < n; k++) {
    TEUCHOS_ASSERT_EQUALITY(0, innerY(k)->Update(lambda[k], rightBordering, 1.0));
  }

  // Add lower bordering.
  std::vector<double> alpha(n);
  TEUCHOS_ASSERT_EQUALITY(0, lowerBordering.Dot(innerX, &alpha[0]));
  for (int k = 0; k < n; k++) {
    alpha[k] += lambda[k] * d_;
  }

  // Merge it all together.
  Nosh::BorderingHelpers::merge(innerY, &alpha[0], Y);

  return 0;
}
// =============================================================================
int
BorderedOperator::
ApplyInverse(const Tpetra::MultiVector<double,int,int> &X,
             Tpetra::MultiVector<double,int,int> &Y
            ) const
{
  // Inverse via Schur formulation.
  //
  // [A     B]^{-1} [X]
  // [<C,.> D]      [lambda]
  // =
  // [A^{-1} X + A^{-1} B S^{-1} <C, A^{-1} X> - A^{-1} B S^{-1} lambda]
  // [-S^{-1} <C, A^{-1} X>                    + S^{-1} lambda         ].
#ifndef NDEBUG
  TEUCHOS_ASSERT(X.getMap().SameAs(domainMap_));
  TEUCHOS_ASSERT(Y.getMap().SameAs(rangeMap_));
#endif
  const int n = X.NumVectors();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(n, Y.NumVectors());
#endif
  // Dissect X.
  Tpetra::Vector<double,int,int> innerX(innerOperator_->getDomainMap());
  std::vector<double> lambda(n);
  Nosh::BorderingHelpers::dissect(X, innerX, &lambda[0]);
  // Apply inverse inner operator with right hand side `right bordering'.
  // TODO useTranspose_
  Tpetra::Vector<double,int,int> AiB(innerOperator_->getRangeMap());
  TEUCHOS_ASSERT_EQUALITY(0, innerOperator_->ApplyInverse(b_, AiB));

  // Schur complement S = D - <C, A^{-1} B>.
  double s;
  TEUCHOS_ASSERT_EQUALITY(0, c_.Dot(AiB, &s));
  s = d_ - s;
  TEUCHOS_ASSERT_INEQUALITY(fabs(s), >=, 1.0e-15);

  Tpetra::MultiVector<double,int,int> innerY(innerOperator_->getRangeMap(), n);
  Tpetra::Vector<double,int,int> AiX(innerOperator_->getRangeMap());
  std::vector<double> alpha(n);
  for (int k = 0; k < n; k++) {
    // innerY = A^{-1} X
    TEUCHOS_ASSERT_EQUALITY(0, innerOperator_->ApplyInverse(innerX, *(innerY(k))));
    double cAiX;
    TEUCHOS_ASSERT_EQUALITY(0, c_.Dot(AiX, &cAiX));

    // [A^{-1} X + A^{-1} B S^{-1} <C, A^{-1} X> - A^{-1} B S^{-1} lambda]
    innerY(k)->Update((cAiX-lambda[k])/s, AiB, 1.0);

    // [-S^{-1} <C, A^{-1} X>                    + S^{-1} lambda         ].
    alpha[k] = (-cAiX + lambda[k]) / s;
  }

  // Merge it all together.
  Nosh::BorderingHelpers::merge(innerY, &alpha[0], Y);

  return 0;
}
// =============================================================================
double
BorderedOperator::
NormInf() const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(false,
                              "Not yet implemented.");
  return 0.0;
}
// =============================================================================
const char *
BorderedOperator::
Label() const
{
  return "Bordered operator";
}
// =============================================================================
bool
BorderedOperator::
UseTranspose() const
{
  return useTranspose_;
}
// =============================================================================
bool
BorderedOperator::
HasNormInf() const
{
  return false;
}
// =============================================================================
const Teuchos::Comm<int> &
BorderedOperator::
Comm() const
{
  return innerOperator_->Comm();
}
// =============================================================================
const Tpetra::Map<int,int> &
BorderedOperator::
getDomainMap() const
{
  return domainMap_;
}
// =============================================================================
const Tpetra::Map<int,int> &
BorderedOperator::
getRangeMap() const
{
  return rangeMap_;
}
// =============================================================================
const std::shared_ptr<Tpetra::Operator<double,int,int>>
BorderedOperator::
getInnerOperator() const
{
  return innerOperator_;
}
// =============================================================================
void
BorderedOperator::
resetBordering(const Tpetra::Vector<double,int,int> & b,
               const Tpetra::Vector<double,int,int> & c,
               const double d
              )
{
  b_ = b;
  c_ = c;
  d_ = d;
  return;
}
// =============================================================================
} // namespace Nosh
