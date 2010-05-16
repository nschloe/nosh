/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl"omer

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

#ifndef GINLA_STATE_H
#define GINLA_STATE_H

#include "Ginla_Typedefs.h"
#include "Recti_Grid_General.h"

namespace Ginla {

class State
{
public:
  
  //! Constructor.
  State( const Teuchos::RCP<ComplexVector>              & psi,
         const Teuchos::RCP<const Recti::Grid::General> & grid );

  //! Constructor without \f$\psi\f$. The values will be initialized to 0.
  State( const Teuchos::RCP<const ComplexMap>           & map,
         const Teuchos::RCP<const Recti::Grid::General> & grid );
  
  //! Const getter.
  Teuchos::RCP<const ComplexVector>
  getValuesConst () const;
  
  //! Nonconst getter for the values.
  Teuchos::RCP<ComplexVector>
  getValuesNonConst ();
  
  double
  getChi () const;
  
  void
  setValues ( Teuchos::RCP<ComplexVector> & psi );
  
  const Teuchos::RCP<const Recti::Grid::General>
  getGrid () const;
  
  //! Save the state to file \c fileName together with the parameters \c p.
  void
  save( const std::string            & fileName,
        const Teuchos::ParameterList & p
      ) const;
  
  //! Just plain save the file to \c fileName.
  void
  save( const std::string & fileName
      ) const;
         
  //! \f$L^2(\Omega)\f$-norm.
  double
  normalizedScaledL2Norm () const;
  
  //! Updates the values of the state according to
  //! \f$a \leftarrow \alpha a + \beta b \f$.
  void
  update( const double alpha,
          const Ginla::State b,
          const double beta  );

  /** Calcuate the grid approximation of the Gibbs free energy
    * \f[
    * \mathcal{G} = \int\nolimits_{\Omega} |\psi|^4 \,\mathrm{d}\omega
    * \f]
    * of a given state \f$\psi\f$.
    */
  double
  freeEnergy () const;

  /*! Calculate the vorticity of the current solution. */
  int
  getVorticity () const;
                 
  /** Writes a solution \c psi to a file with all parameters that
    * may be interesting.
    */
  void
  writeStateToFile ( LOCA::ParameterVector & params,
                     const std::string     & fileBaseName,
                     const std::string     & outputFormat
                   );

protected:
private:
  
  //! Numerical values.
  Teuchos::RCP<ComplexVector> psi_;
  

  double chi_;
  
  //! The grid on which the state exists.
  const Teuchos::RCP<const Recti::Grid::General> grid_;
  
};

// nonmember functions

}

#endif // GINLA_STATE_H
