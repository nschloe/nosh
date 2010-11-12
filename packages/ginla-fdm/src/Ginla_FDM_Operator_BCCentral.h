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

#ifndef GINLA_FDM_OPERATOR_BCCENTRAL_H
#define GINLA_FDM_OPERATOR_BCCENTRAL_H

#include "Ginla_FDM_Operator_Virtual.h"

namespace Ginla {
  namespace FDM {
  namespace Operator {

class BCCentral:
            public Virtual
{
public:

    //! Default constructor.
    BCCentral ( const Teuchos::RCP<Recti::Grid::Uniform>                    & grid,
                const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & A,
                const Teuchos::RCP<const ComplexMap>                        & domainMap,
                const Teuchos::RCP<const ComplexMap>                        & rangeMap
              );

    //! Destructor
    virtual
    ~BCCentral();

    virtual
    Teuchos::RCP<Ginla::State::Updatable>
    getF( const Teuchos::RCP<const Ginla::FDM::State> & state
        ) const;

    virtual
    Teuchos::RCP<const Ginla::State::Updatable>
    getDFDh0( const Teuchos::RCP<const Ginla::FDM::State> & state
            ) const;

    virtual Teuchos::RCP<const Komplex2::DoubleMatrix>
    getJacobian ( const Teuchos::RCP<const Ginla::FDM::State> & state );

protected:
private:
    //! Return the value of the Ginzburg-Landau equations for the equation
    //! at eqType.
    double_complex
    getFEntry_ ( Teuchos::ArrayRCP<const double_complex> & psiView,
                 const double                              chi,
                 const int localIndex,
                 const int globalIndex,
                 Teuchos::ArrayRCP<const double> & ALeftView,
                 Teuchos::ArrayRCP<const double> & ARightView,
                 Teuchos::ArrayRCP<const double> & ABelowView,
                 Teuchos::ArrayRCP<const double> & AAboveView,
                 const double                      h
               ) const;

    //! Return the value of the Ginzburg-Landau equations for the equation
    //! at eqType.
    double_complex
    getDFDh0Entry_( Teuchos::ArrayRCP<const double_complex> & psiView,
                    const double                              chi,
                    const int localIndex,
                    const int globalIndex,
                    Teuchos::ArrayRCP<const double> & ALeftView,
                    Teuchos::ArrayRCP<const double> & ARightView,
                    Teuchos::ArrayRCP<const double> & ABelowView,
                    Teuchos::ArrayRCP<const double> & AAboveView,
                    Teuchos::ArrayRCP<const double> & dAdH0LeftView,
                    Teuchos::ArrayRCP<const double> & dAdH0RightView,
                    Teuchos::ArrayRCP<const double> & dAdH0BelowView,
                    Teuchos::ArrayRCP<const double> & dAdH0AboveView,
                    const double                      h
                  ) const;

    //! Returns entries and positions of the Jacobian matrix belonging to the
    //! boundary conditions.
    void
    getJacobianRow_ ( Teuchos::ArrayRCP<const double_complex> & psiView,
                      const double                              chi,
                      const int localIndex,
                      const int globalIndex,
                      Teuchos::ArrayRCP<const double> & ALeftView,
                      Teuchos::ArrayRCP<const double> & ARightView,
                      Teuchos::ArrayRCP<const double> & ABelowView,
                      Teuchos::ArrayRCP<const double> & AAboveView,
                      const double                      h,
                      Teuchos::Array<Thyra::Ordinal>  & columnIndicesPsi,
                      Teuchos::Array<double_complex>  & valuesPsi,
                      Teuchos::Array<Thyra::Ordinal>  & columnIndicesPsiConj,
                      Teuchos::Array<double_complex>  & valuesPsiConj
                    ) const;

};

    } // namespace Operator
  }
} // namespace GL

#endif // GINLA_FDM_OPERATOR_BCCENTRAL_H
