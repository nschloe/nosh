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

#ifndef NOSH_KEOPRECONDITIONER_H
#define NOSH_KEOPRECONDITIONER_H
// =============================================================================
#include <Epetra_Vector.h>
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
class StkMesh;
namespace ScalarField {
class Virtual;
}
namespace MatrixBuilder {
class Virtual;
}
}
namespace Belos {
class EpetraPrecOp;
}
class Amesos_BaseSolver;
class Epetra_LinearProblem;
namespace ML_Epetra {
class MultiLevelPreconditioner;
}
// =============================================================================
namespace Nosh {
// =============================================================================
class KeoRegularized : public Epetra_Operator
{
public:
KeoRegularized(const Teuchos::RCP<const Nosh::StkMesh> &mesh,
               const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness,
               const Teuchos::RCP<const Nosh::MatrixBuilder::Virtual> &matrixBuilder
               );

// Destructor.
~KeoRegularized();

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

void
rebuild(const double g,
        const Teuchos::Array<double> &mvpParams,
        const Teuchos::RCP<const Epetra_Vector> &psi
       );

protected:
private:

void
rebuildInverse_();

void
rebuildAbsPsiSquared_(const Teuchos::RCP<const Epetra_Vector> &psi);

private:

private:
bool useTranspose_;

const Teuchos::RCP<const Nosh::StkMesh> mesh_;
double g_;
const Teuchos::RCP<const Nosh::ScalarField::Virtual> thickness_;
const Teuchos::RCP<const Nosh::MatrixBuilder::Virtual> matrixBuilder_;

// |psi|^2
Epetra_Vector absPsiSquared_;

// Make sure to create the matrix in memory only once and then
// override it as necessary. The reason for this is that ML
// gets initialized only once and, upon ML.recompute(), relies
// on the (new) data being available at the same adress.
// Failure to comply to this will lead to memory errors.
Epetra_FECrsMatrix regularizedMatrix_;

const Epetra_Comm &comm_;

Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MlPrec_;
const int numCycles_;

#ifdef NOSH_TEUCHOS_TIME_MONITOR
const Teuchos::RCP<Teuchos::Time> timerRebuild0_;
const Teuchos::RCP<Teuchos::Time> timerRebuild1_;
#endif

Teuchos::RCP<Teuchos::FancyOStream> out_;
};
// =============================================================================
} // namespace Nosh

#endif // NOSH_KEOPRECONDITIONER_H