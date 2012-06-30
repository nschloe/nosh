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

#ifndef GINLA_KEOPRECONDITIONER_H
#define GINLA_KEOPRECONDITIONER_H
// =============================================================================
// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include "Ginla_config.h"

#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  #include <Teuchos_Time.hpp>
#endif
#include <Teuchos_FancyOStream.hpp>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
// =============================================================================
namespace Ginla {
class KeoContainer;
class StkMesh;
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
namespace Ginla {
// =============================================================================
class KeoRegularized : public Epetra_Operator
{
public:
KeoRegularized( const Teuchos::RCP<const Ginla::StkMesh> &mesh,
                const double g,
                const Teuchos::RCP<const Epetra_Vector> &thickness,
                const Teuchos::RCP<Ginla::KeoContainer> &keoContainer,
                const Teuchos::RCP<const Epetra_Vector> &psi
              );

// Destructor.
~KeoRegularized();

virtual int
SetUseTranspose( bool UseTranspose );

virtual int
Apply( const Epetra_MultiVector &X,
       Epetra_MultiVector &Y
       ) const;

virtual int
ApplyInverse( const Epetra_MultiVector &X,
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
rebuild(const Teuchos::Array<double> &mvpParams,
        const Teuchos::RCP<const Epetra_Vector> &psi
       );

void
rebuildInverse();

protected:
private:

void
rebuildAbsPsiSquared_(const Teuchos::RCP<const Epetra_Vector> &psi);

private:

private:
bool useTranspose_;

const Teuchos::RCP<const Ginla::StkMesh> mesh_;
const double g_;
const Teuchos::RCP<const Epetra_Vector> thickness_;
Teuchos::RCP<Ginla::KeoContainer> keoContainer_;

// |psi|^2
const Teuchos::RCP<Epetra_Vector> absPsiSquared_;

// Make sure the matrix is persistent in memory. ML requires that.
Teuchos::RCP<Epetra_CrsMatrix> keoRegularizedMatrix_;

const Epetra_Comm &comm_;

Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MlPrec_;
const int numCycles_;

#ifdef GINLA_TEUCHOS_TIME_MONITOR
const Teuchos::RCP<Teuchos::Time> timerRebuild0_;
const Teuchos::RCP<Teuchos::Time> timerRebuild1_;
#endif

Teuchos::RCP<Teuchos::FancyOStream> out_;
};
// =============================================================================
} // namespace Ginla

#endif // GINLA_KEOPRECONDITIONER_H
