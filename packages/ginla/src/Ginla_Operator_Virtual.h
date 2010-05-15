/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl"omer

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

#include "Ginla_Typedefs.h"
#include "Ginla_State.h"

// forward declarations
namespace Recti {
  namespace Grid {
    class Uniform;
    class General;
  }
}
namespace Ginla {
  namespace MagneticVectorPotential{
    class Centered;
  }
}

#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

namespace Ginla {
  
  namespace Operator {

class Virtual
{
public:
    //! Default constructor.
    Virtual ( const Teuchos::RCP<Recti::Grid::Uniform>                     & grid,
              const Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> & A,
              const Teuchos::RCP<const ComplexMap>                         & domainMap,
              const Teuchos::RCP<const ComplexMap>                         & rangeMap
            );

    //! Destructor
    virtual ~Virtual();
    
    virtual Teuchos::RCP<Ginla::State>
    getF( const Teuchos::RCP<const Ginla::State> & state ) const = 0; // purely virtual

    virtual void
    getJacobianRow ( const Teuchos::RCP<const Ginla::State> & state,
                     const int                                k,
                     Teuchos::Array<int>                    & columnIndicesPsi,
                     Teuchos::Array<double_complex>         & valuesPsi,
                     Teuchos::Array<int>                    & columnIndicesPsiConj,
                     Teuchos::Array<double_complex>         & valuesPsiCon
                   ) const = 0; // purely virtual

    Teuchos::RCP<const Recti::Grid::General>
    getGrid() const;

    void
    setParameters( const LOCA::ParameterVector & p );
    
    Teuchos::RCP<LOCA::ParameterVector>
    getParameters() const;

protected:
    const Teuchos::RCP<const ComplexMap> domainMap_;
    const Teuchos::RCP<const ComplexMap> rangeMap_;
    const Teuchos::RCP<Recti::Grid::Uniform> grid_;
    const Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> A_;

private:
};

  } // namespace Virtual
  
} // namespace GL
#endif // GLOPERATORVIRTUAL_H
