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

#ifndef GINLA_FVM_KINETICENERGYOPERATOR_H
#define GINLA_FVM_KINETICENERGYOPERATOR_H
// =============================================================================
#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>

#include "Ginla_MagneticVectorPotential_Virtual.h"
// =============================================================================
namespace Ginla {
namespace FVM {
// =============================================================================
class KineticEnergyOperator : public Epetra_Operator
{
public:
    KineticEnergyOperator();

    // Destructor.
    ~KineticEnergyOperator();

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

    void
    setParameters( const double mu,
                   const Teuchos::Tuple<double,3> & scaling
                 );

protected:
private:
    const Teuchos::RCP<VIO::Mesh::Mesh> mesh_;

    const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp_;

    Teuchos::RCP<Epetra_FECrsGraph> keoGraph_;
    Teuchos::RCP<Epetra_FECrsMatrix> keo_;

    double mu_;
    Teuchos::Tuple<double,3> scaling_;

    double keoMu_;
    Teuchos::Tuple<double,3> keoScaling_;
};
// =============================================================================
} // namespace FVM
} // namespace Ginla

#endif // GINLA_FVM_KINETICENERGYOPERATOR_H
