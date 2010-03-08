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

#include "GL_MagneticVectorPotential_Centered.h"
#include "Recti_Grid_Uniform.h"

typedef std::complex<double> double_complex;

namespace GL {
  
  namespace Operator {

class Virtual
{
public:
    //! Default constructor.
    Virtual ( Teuchos::RCP<Recti::Grid::Uniform>                  & grid,
              Teuchos::RCP<GL::MagneticVectorPotential::Centered> & A
            );

    //! Destructor
    virtual ~Virtual();

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

    Teuchos::RCP<const Recti::Grid::General>
    getGrid() const;

    void
    setParameters( const LOCA::ParameterVector & p );
    
    Teuchos::RCP<LOCA::ParameterVector>
    getParameters() const;

protected:
    Teuchos::RCP<const ComplexVector> psi_;
    double chi_;
    Teuchos::RCP<Recti::Grid::Uniform> grid_;
    Teuchos::RCP<GL::MagneticVectorPotential::Centered> A_;

private:
};

  } // namespace Virtual
  
} // namespace GL
#endif // GLOPERATORVIRTUAL_H
