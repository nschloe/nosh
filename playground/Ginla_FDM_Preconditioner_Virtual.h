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

#ifndef GINLA_PRECONDITIONER_VIRTUAL_H
#define GINLA_PRECONDITIONER_VIRTUAL_H

#include "Komplex2_DoubleMatrix.h"
#include "Ginla_State.h"

// forward declarations
namespace Recti {
  namespace Grid {
    class Uniform;
  }
}

namespace Ginla {

namespace Preconditioner {

class Virtual
{
public:
    //! Default constructor.
    Virtual (  const Teuchos::RCP<Recti::Grid::Uniform> & grid,
               const Teuchos::RCP<const ComplexMap>     & domainMap,
               const Teuchos::RCP<const ComplexMap>     & rangeMap
             );

    //! Destructor
    virtual ~Virtual();
    
    virtual Teuchos::RCP<const Komplex2::DoubleMatrix>
    getMatrix ( const Teuchos::RCP<const Ginla::State> & state ) = 0; // purely virtual

protected:
  
  const Teuchos::RCP<const ComplexMap> domainMap_;
  const Teuchos::RCP<const ComplexMap> rangeMap_;
  const Teuchos::RCP<Recti::Grid::Uniform> grid_;
  
  //! Data for building the double matrix A/B, corresponding to \f$\psi\f$ and \f$\overline{\psi}\f$, respectively.
  bool firstTime_;
  Teuchos::RCP<Komplex2::DoubleMatrix> AB_;
    
private:
  
};

}

}

#endif // GINLA_PRECONDITIONER_VIRTUAL_H
