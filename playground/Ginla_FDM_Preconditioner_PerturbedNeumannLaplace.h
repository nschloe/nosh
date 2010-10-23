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

#ifndef GINLA_PRECONDITIONER_PERTURBEDNEUMANNLAPLACE_H
#define GINLA_PRECONDITIONER_PERTURBEDNEUMANNLAPLACE_H

#include "Ginla_Preconditioner_Virtual.h"
#include "Ginla_Typedefs.h"

#include <Teuchos_RCP.hpp>

// forward declarations
namespace Ginla {
  class State;
  namespace Komplex {
    class DoubleMatrix;
  }
}

namespace Recti {
  namespace Grid {
    class Uniform;
  }
}

namespace Ginla {

namespace Preconditioner {

class PerturbedNeumannLaplace: public Ginla::Preconditioner::Virtual
{
public:
    //! Default constructor.
    PerturbedNeumannLaplace ( const Teuchos::RCP<Recti::Grid::Uniform> & grid,
                              const Teuchos::RCP<const ComplexMap>     & domainMap,
                              const Teuchos::RCP<const ComplexMap>     & rangeMap
                            );

    //! Destructor
    virtual ~PerturbedNeumannLaplace();
    
    Teuchos::RCP<const Komplex2::DoubleMatrix>
    getMatrix ( const Teuchos::RCP<const Ginla::State> & state );
    
protected:
  
private:
  const double epsilon;
  
private:

  void
  getMatrixRow_ ( const int localIndex,
                  const int globalIndex,
                  const double                      h,
                  Teuchos::Array<Thyra::Ordinal>  & columnIndicesPsi,
                  Teuchos::Array<double_complex>  & valuesPsi,
                  Teuchos::Array<Thyra::Ordinal>  & columnIndicesPsiConj,
                  Teuchos::Array<double_complex>  & valuesPsiConj
                ) const;

};

}

}

#endif // GINLA_PRECONDITIONER_PERTURBEDNEUMANNLAPLACE_H
