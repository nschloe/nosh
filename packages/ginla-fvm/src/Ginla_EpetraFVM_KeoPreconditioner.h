/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef GINLA_EPETRAFVM_KEOPRECONDITIONER_H
#define GINLA_EPETRAFVM_KEOPRECONDITIONER_H
// =============================================================================
// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Time.hpp>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_LinearProblem.h>

#include "Ginla_MagneticVectorPotential_Virtual.h"
// =============================================================================
namespace Ginla {
    namespace EpetraFVM {
        class KeoFactory;
        class StkMesh;
    }
}
namespace Belos {
    class EpetraPrecOp;
}
class Amesos_BaseSolver;
class Epetra_LinearProblem;
// =============================================================================
namespace Ginla {
namespace EpetraFVM {
// =============================================================================
class KeoPreconditioner: public Epetra_Operator
{
public:
    KeoPreconditioner( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh>               & mesh,
                       const Teuchos::RCP<const Epetra_Vector>                     & thickness,
                       const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp
                     );

    // Destructor.
    ~KeoPreconditioner();

    virtual int
    SetUseTranspose( bool UseTranspose );

    virtual int
    Apply ( const Epetra_MultiVector & X,
                  Epetra_MultiVector & Y
          ) const;

    virtual int
    ApplyInverse ( const Epetra_MultiVector & X,
                         Epetra_MultiVector & Y
                 ) const;

    int
    ApplyInverse2 ( const Epetra_MultiVector & X,
                         Epetra_MultiVector & Y
                 ) const;

    virtual double
    NormInf () const;

    virtual const char *
    Label () const;

    virtual bool
    UseTranspose () const;

    virtual bool
    HasNormInf () const;

    virtual const Epetra_Comm &
    Comm () const;

    virtual const Epetra_Map & 	OperatorDomainMap () const;

    virtual const Epetra_Map & 	OperatorRangeMap () const;

public:

    void
    rebuild();

    void
    rebuild( const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams,
             const Teuchos::Tuple<double,3>                  & scaling
           );

protected:
private:
    int
    ApplyInverseMl_( const Epetra_MultiVector & X,
                           Epetra_MultiVector & Y
                   ) const;

    int
    ApplyInverseIlu_ ( const Epetra_MultiVector & X,
                             Epetra_MultiVector & Y
                     ) const;

    void
    rebuildMl_();

    void
    rebuildIlu_();

private:

    enum InversionType { Ilu, Ml };

private:
    bool useTranspose_;
    const Teuchos::RCP<const Epetra_Comm> comm_;

    Teuchos::RCP<Ginla::EpetraFVM::KeoFactory> keoFactory_;
    Teuchos::RCP<Epetra_FECrsMatrix> keoPrec_;

    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec_;

    Teuchos::RCP<Epetra_LinearProblem> keoProblem_;
    Teuchos::RCP<Amesos_BaseSolver> keoSolver_;

    InversionType invType_;

    const Teuchos::RCP<Teuchos::Time> rebuildTime_;
    const Teuchos::RCP<Teuchos::Time> rebuildMlTime_;
    const Teuchos::RCP<Teuchos::Time> rebuildIluTime_;
    const Teuchos::RCP<Teuchos::Time> regularizationTime_;
};
// =============================================================================
} // namespace FVM
} // namespace Ginla

#endif // GINLA_EPETRAFVM_KEOPRECONDITIONER_H
