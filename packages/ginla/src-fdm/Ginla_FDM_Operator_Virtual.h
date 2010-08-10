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

#ifndef GINLA_FDM_OPERATOR_VIRTUAL_H
#define GINLA_FDM_OPERATOR_VIRTUAL_H

#include "Ginla_Typedefs.h"
#include "Ginla_FDM_State.h"
#include "Ginla_ParameterHost_Virtual.h"

// forward declarations
namespace Recti {
  namespace Grid {
    class Uniform;
    class General;
  }
}
namespace Ginla {
  namespace MagneticVectorPotential{
    class Virtual;
  }
  namespace Komplex {
    class DoubleMatrix; 
  }
}

#include <Teuchos_Array.hpp>
#include <Epetra_Vector.h>
#include <LOCA_Parameter_Vector.H>

namespace Ginla {
  namespace FDM {
  namespace Operator {

class Virtual:
    public Ginla::ParameterHost::Virtual
{
public:
    //! Default constructor.
    Virtual ( const Teuchos::RCP<Recti::Grid::Uniform>                    & grid,
              const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & A,
              const Teuchos::RCP<const ComplexMap>                        & domainMap,
              const Teuchos::RCP<const ComplexMap>                        & rangeMap
            );

    //! Destructor
    virtual ~Virtual();
    
    virtual Teuchos::RCP<Ginla::State::Virtual>
    getF( const Teuchos::RCP<const Ginla::FDM::State> & state ) const = 0; // purely virtual
    
    virtual Teuchos::RCP<const Ginla::State::Virtual>
    getDFDh0( const Teuchos::RCP<const Ginla::FDM::State> & state ) const = 0; // purely virtual

    virtual Teuchos::RCP<const Ginla::Komplex::DoubleMatrix>
    getJacobian ( const Teuchos::RCP<const Ginla::FDM::State> & state ) = 0; // purely virtual

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
    
    //! Builds the cache for values of \f$A\f$ and \f$\frac{\text{d}A}{\text{d}H0}\f$.
    void
    buildACache_() const;
    
    mutable bool cacheNeedsUpdating_;
    
    // cache for the queries to A
    const Teuchos::RCP<RealVector> ALeft_;
    const Teuchos::RCP<RealVector> ARight_;
    const Teuchos::RCP<RealVector> AAbove_;
    const Teuchos::RCP<RealVector> ABelow_;
    
    // cache for the queries to dAdH0
    const Teuchos::RCP<RealVector> dAdH0Left_;
    const Teuchos::RCP<RealVector> dAdH0Right_;
    const Teuchos::RCP<RealVector> dAdH0Above_;
    const Teuchos::RCP<RealVector> dAdH0Below_;

    //! Data for building the double matrix A/B, corresponding to \f$\psi\f$ and \f$\overline{\psi}\f$, respectively.
    bool firstTime_;
    Teuchos::RCP<Ginla::Komplex::DoubleMatrix> AB_;
    
private:
  
    const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> A_;

};

    } // namespace Virtual
  }
} // namespace GL
#endif // GINLA_FDM_OPERATOR_VIRTUAL_H
