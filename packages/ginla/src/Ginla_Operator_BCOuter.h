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

#ifndef GLOPERATORBCOUTER_H
#define GLOPERATORBCOUTER_H

#include "Ginla_Operator_Virtual.h"

namespace Ginla {
  namespace Operator {

class BCOuter:
            public Virtual
{
public:

    //! Default constructor.
    BCOuter ( const Teuchos::RCP<Recti::Grid::Uniform>             & grid,
              const Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> & A
            );

    //! Destructor
    virtual
    ~BCOuter();

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

  } // namespace Operator
} // namespace GL

#endif // GLOPERATORBCOUTER_H
