/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

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

#ifndef GLOPERATORVIRTUAL_H
#define GLOPERATORVIRTUAL_H

#include <complex>
#include <Tpetra_Vector.hpp>

#include "MagneticVectorPotential.h"
#include "GridUniform.h"

typedef std::complex<double> double_complex;

class GlOperatorVirtual
{
public:
    //! Default constructor.
    GlOperatorVirtual ( Teuchos::RCP<GridUniform>             & grid,
                        Teuchos::RCP<MagneticVectorPotential> & A
                      );

    //! Destructor
    virtual ~GlOperatorVirtual();

    virtual double_complex
    getEntry ( const int k ) const = 0; // purely virtual

    virtual void
    getJacobianRow ( const int                        k,
                     Teuchos::Array<int>            & columnIndicesPsi,
                     Teuchos::Array<double_complex> & valuesPsi,
                     Teuchos::Array<int>            & columnIndicesPsiConj,
                     Teuchos::Array<double_complex> & valuesPsiCon
                   ) const = 0; // purely virtual

    void
    updatePsi ( const Teuchos::RCP<const ComplexVector> psi );

    void
    setChi ( const double chi );

    void
    setH0 ( const double h0 );

    double
    getH0 () const;

    Teuchos::RCP<const Grid>
    getGrid() const;

    double
    getScaling () const;

    void
    setScaling ( const double scaling );

protected:
    Teuchos::RCP<const ComplexVector> psi_;
    double chi_;
    Teuchos::RCP<GridUniform> grid_;
    Teuchos::RCP<MagneticVectorPotential> A_;

private:
};

#endif // GLOPERATORVIRTUAL_H
