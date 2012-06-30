// @HEADER
//
//    Ginla Jacobian operator.
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

#ifndef GINLA_JACOBIANOPERATOR_H
#define GINLA_JACOBIANOPERATOR_H
// =============================================================================
// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>
// =============================================================================
// forward declarations
namespace Ginla {
class KeoContainer;
class StkMesh;
namespace ScalarPotential{
class Virtual;
}
}
class Epetra_CrsGraph;
class Epetra_FECrsMatrix;
class Epetra_Vector;
// =============================================================================
namespace Ginla {

class JacobianOperator : public Epetra_Operator
{
public:
JacobianOperator(const Teuchos::RCP<const Ginla::StkMesh> &mesh,
                 const Teuchos::RCP<Ginla::ScalarPotential::Virtual> &scalarPotential,
                 const double g,
                 const Teuchos::RCP<const Epetra_Vector> &thickness,
                 const Teuchos::RCP<Ginla::KeoContainer> &keoContainer,
                 const Teuchos::RCP<const Epetra_Vector> &current_X
                 );

// Destructor.
~JacobianOperator ();

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
rebuild( const Teuchos::Array<double> &mvpParams,
         const Teuchos::Array<double> &spParams,
         const Teuchos::RCP<const Epetra_Vector> &current_X
         );

protected:

private:
void
rebuildDiags_(const Teuchos::Array<double> &spParams,
              const Teuchos::RCP<const Epetra_Vector> &current_X
              );

private:
bool useTranspose_;

const Teuchos::RCP<const Ginla::StkMesh> mesh_;
const Teuchos::RCP<Ginla::ScalarPotential::Virtual> scalarPotential_;
const double g_;
const Teuchos::RCP<const Epetra_Vector> thickness_;
const Teuchos::RCP<Ginla::KeoContainer> keoContainer_;
Teuchos::RCP<const Epetra_FECrsMatrix> keo_;

Teuchos::RCP<const Epetra_Vector> current_X_;

const Teuchos::RCP<Epetra_Vector> diag0_;
const Teuchos::RCP<Epetra_Vector> diag1b_;
};

} // namespace Ginla
// =============================================================================
#endif // GINLA_JACOBIANOPERATOR_H
