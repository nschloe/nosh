// @HEADER
//
//    Regularized kinetic energy operator.
//    Copyright (C) 2010--2012  Nico Schlömer
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

#ifndef NOSH_BORDEREDOPERATOR_H
#define NOSH_BORDEREDOPERATOR_H
// =============================================================================
#include <Tpetra_Vector.hpp>
#include <Tpetra::Map<int,int>.h>
#include <Tpetra::Operator<double,int,int>.h>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Teuchos_FancyOStream.hpp>
// =============================================================================
namespace nosh
{
// =============================================================================
class BorderedOperator : public Tpetra::Operator<double,int,int>
{
public:
  BorderedOperator(
      const std::shared_ptr<Tpetra::Operator<double,int,int>> & innerOperator,
      const Tpetra::Vector<double,int,int> & b,
      const Tpetra::Vector<double,int,int> & c,
      const double d
      );

  // Destructor.
  ~BorderedOperator();

  virtual int
  SetUseTranspose(bool UseTranspose);

  virtual int
  Apply(const Tpetra::MultiVector<double,int,int> &X,
        Tpetra::MultiVector<double,int,int> &Y
      ) const;

  virtual int
  ApplyInverse(const Tpetra::MultiVector<double,int,int> &X,
               Tpetra::MultiVector<double,int,int> &Y
             ) const;

  virtual double
  NormInf() const;

  virtual const char *
  Label() const;

  virtual bool
  UseTranspose() const;

  virtual bool
  HasNormInf() const;

  virtual const Teuchos::Comm<int> &
  Comm() const;

  virtual const Tpetra::Map<int,int> &getDomainMap() const;

  virtual const Tpetra::Map<int,int> &getRangeMap() const;

public:
  const std::shared_ptr<Tpetra::Operator<double,int,int>>
  getInnerOperator() const;

  void
  resetBordering(const Tpetra::Vector<double,int,int> & b,
                 const Tpetra::Vector<double,int,int> & c,
                 const double d
               );

protected:
private:
  const std::shared_ptr<Tpetra::Operator<double,int,int>> innerOperator_;
  Tpetra::Vector<double,int,int> b_;
  Tpetra::Vector<double,int,int> c_;
  double d_;
  bool useTranspose_;
  const Tpetra::Map<int,int> domainMap_;
  const Tpetra::Map<int,int> rangeMap_;
};
// =============================================================================
} // namespace nosh
#endif // NOSH_BORDEREDOPERATOR_H
