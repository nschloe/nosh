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

#ifndef GLOPERATORBCCENTRAL_H
#define GLOPERATORBCCENTRAL_H

#include "GlOperatorVirtual.h"

class GlOperatorBCCentral:
            public GlOperatorVirtual
{
public:

    //! Default constructor.
    GlOperatorBCCentral ( Teuchos::RCP<GridUniform>             & grid,
                          Teuchos::RCP<MagneticVectorPotential> & A
                        );

    //! Destructor
    virtual
    ~GlOperatorBCCentral();

    //! Return the value of the Ginzburg-Landau equations for the equation
    //! at eqType.
    virtual double_complex
    getEntry ( const int k ) const;

    //! Returns entries and positions of the Jacobian matrix belonging to the
    //! boundary conditions.
    virtual void
    getJacobianRow ( const int                        k,
                     Teuchos::Array<int>            & columnIndicesPsi,
                     Teuchos::Array<double_complex> & valuesPsi,
                     Teuchos::Array<int>            & columnIndicesPsiConj,
                     Teuchos::Array<double_complex> & valuesPsiCon
                   ) const;

protected:
private:

};

#endif // GLOPERATORBCCENTRAL_H
