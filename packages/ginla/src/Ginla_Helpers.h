/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl\"omer

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

#ifndef GL_HELPERS_H
#define GL_HELPERS_H

#include <LOCA_Parameter_Vector.H>
#include <Teuchos_ParameterList.hpp>

#include "Recti_Grid_General.h"

namespace Ginla {

namespace Helpers
{
  //! Method that reads a given Teuchos::ParameterList and puts all \c double
  //! entries into a LOCA::ParameterVector.
  Teuchos::RCP<LOCA::ParameterVector>
  teuchosParameterList2locaParameterVector( const Teuchos::ParameterList & p
                                          );

  Teuchos::RCP<Teuchos::ParameterList>
  locaParameterVector2teuchosParameterList( const LOCA::ParameterVector & pL );
                                          
  //! Merges two \c LOCA::ParameterLists into one, checking for discrepancies
  //! in the entries.
  Teuchos::RCP<LOCA::ParameterVector>
  mergeLocaParameterVectors( const LOCA::ParameterVector & p0,
                             const LOCA::ParameterVector & p1
                           );
                           
  void
  appendToTeuchosParameterList( Teuchos::ParameterList      & p,
                                const LOCA::ParameterVector & pL,
                                const std::string           & labelPrepend = ""
                              );

  double
  normalizedScaledL2Norm ( const ComplexVector        & psi,
                           const Recti::Grid::General & grid );

  /** Calcuate the grid approximation of the Gibbs free energy
    * \f[
    * \mathcal{G} = \int\nolimits_{\Omega} |\psi|^4 \,\mathrm{d}\omega
    * \f]
    * of a given state \f$\psi\f$.
    */
  double
  freeEnergy ( const ComplexVector & psi,
               const Recti::Grid::General & grid );

  /*! Calculate the vorticity of the current solution. */
  int
  getVorticity ( const ComplexVector        & psi,
                 const Recti::Grid::General & grid );
                 
  /** Writes a solution \c psi to a file with all parameters that
    * may be interesting.
    */
  void
  writeStateToFile ( const Teuchos::RCP<const ComplexVector>        & psi,
                     const Teuchos::RCP<const Recti::Grid::General> & grid,
                     LOCA::ParameterVector                          & params,
                     const std::string                              & fileBaseName,
                     const std::string                              & outputFormat
                   );

};

}

#endif // GL_HELPERS_H
