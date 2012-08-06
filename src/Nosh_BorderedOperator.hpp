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

#ifndef NOSH_BORDEREDOPERATOR_H
#define NOSH_BORDEREDOPERATOR_H
// =============================================================================
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  #include <Teuchos_Time.hpp>
#endif
#include <Teuchos_FancyOStream.hpp>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_LinearProblem.h>
// =============================================================================
namespace Nosh {
// =============================================================================
class BorderedOperator : public Epetra_Operator
{
public:
BorderedOperator(const Teuchos::RCP<Epetra_Operator> & innerOperator,
                 const Teuchos::RCP<Epetra_Vector> & b,
                 const Teuchos::RCP<Epetra_Vector> & c,
                 const double d
                 );

// Destructor.
~BorderedOperator();

virtual int
SetUseTranspose( bool UseTranspose );

virtual int
Apply(const Epetra_MultiVector &X,
      Epetra_MultiVector &Y
      ) const;

virtual int
ApplyInverse(const Epetra_MultiVector &X,
             Epetra_MultiVector &Y
             ) const;

virtual double
NormInf() const;

virtual const char *
Label() const;

virtual bool
UseTranspose() const;

virtual bool
HasNormInf() const;

virtual const Epetra_Comm &
Comm() const;

virtual const Epetra_Map &OperatorDomainMap() const;

virtual const Epetra_Map &OperatorRangeMap() const;

public:

const Teuchos::RCP<Epetra_Operator>
getInnerOperator() const;

protected:
private:

const Teuchos::RCP<Epetra_Operator> innerOperator_;
const Teuchos::RCP<Epetra_Vector> b_;
const Teuchos::RCP<Epetra_Vector> c_;
const double d_;
bool useTranspose_;
const Epetra_Map domainMap_;
const Epetra_Map rangeMap_;
};
// =============================================================================
} // namespace Nosh

#endif // NOSH_BORDEREDOPERATOR_H
