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
#include <Teuchos_Tuple.hpp>
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  #include <Teuchos_Time.hpp>
#endif
#include <Teuchos_FancyOStream.hpp>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <LOCA_Parameter_Vector.H>
// =============================================================================
namespace Ginla {
class KeoFactory;
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
                const Teuchos::RCP<const Epetra_Vector> &thickness,
                const Teuchos::RCP<Ginla::KeoFactory> &keoFactory
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
rebuild(const Teuchos::RCP<const Epetra_Vector> &psi);

void
rebuild(const Teuchos::RCP<const LOCA::ParameterVector> &mvpParams,
        const Teuchos::RCP<const Epetra_Vector> &psi
        );

protected:
private:
int
ApplyInverseMl_( const Epetra_MultiVector &X,
                 Epetra_MultiVector &Y
                 ) const;

int
ApplyInverseIlu_( const Epetra_MultiVector &X,
                  Epetra_MultiVector &Y
                  ) const;

void
rebuildMl_();

void
rebuildIlu_();

private:

enum EInversionType { INVERT_ILU, INVERT_ML };

private:
bool useTranspose_;
const Teuchos::RCP<const Ginla::StkMesh> mesh_;
const Teuchos::RCP<const Epetra_Vector> thickness_;

const Epetra_Comm &comm_;

Teuchos::RCP<Ginla::KeoFactory> keoFactory_;
// Make sure the matrix pointer is never changed; ML's
// preconditioner generation depends on that.
//const Teuchos::RCP<Epetra_CrsMatrix> keoRegularized_;
Teuchos::RCP<Epetra_CrsMatrix> keoRegularized_;

Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MlPrec_;

Teuchos::RCP<Epetra_LinearProblem> keoIluProblem_;
Teuchos::RCP<Amesos_BaseSolver> keoIluSolver_;

EInversionType invType_;

#ifdef GINLA_TEUCHOS_TIME_MONITOR
const Teuchos::RCP<Teuchos::Time> timerRebuild_;
#endif

Teuchos::RCP<Teuchos::FancyOStream> out_;
};
// =============================================================================
} // namespace Ginla

#endif // GINLA_KEOPRECONDITIONER_H
