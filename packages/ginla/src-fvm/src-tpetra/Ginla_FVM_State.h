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

#ifndef GINLA_FVM_STATE_H
#define GINLA_FVM_STATE_H

#include "VIO_TpetraMesh_Mesh.h"
#include "Ginla_State_Virtual.h"

namespace LOCA {
  class ParameterVector;
}

namespace Ginla {
namespace FVM {
class State:
  public Ginla::State::Virtual
{
public:

  //! Constructor.
  State( const Teuchos::RCP<ComplexMultiVector>    & psi,
         const Teuchos::RCP<const VIO::TpetraMesh::Mesh> & mesh
       );

  //! Constructor.
  State( const Teuchos::RCP<ComplexVector>         & psi,
         const Teuchos::RCP<const VIO::TpetraMesh::Mesh> & mesh
       );

  //! Constructor without \f$\psi\f$. The values will be initialized to 0.
  State( const Teuchos::RCP<const ComplexMap>      & map,
         const Teuchos::RCP<const VIO::TpetraMesh::Mesh> & mesh
       );

  //! Constructor solely with comminicator and grid. The values will be initialized to 0.
  State( const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
         const Teuchos::RCP<const VIO::TpetraMesh::Mesh>     & mesh
       );

  //! Const getter.
  Teuchos::RCP<const ComplexVector>
  getPsi () const;

  //! Nonconst getter for the values.
  Teuchos::RCP<ComplexVector>
  getPsiNonConst ();

  double
  getChi () const;

  const Teuchos::RCP<const VIO::TpetraMesh::Mesh>
  getMesh () const;

  //! Save the state to file \c fileName together with the parameters \c p.
  void
  save( const std::string            & fileName,
        const Teuchos::ParameterList & p
      ) const;

  //! Just plain save the file to \c fileName.
  void
  save( const std::string & fileName
      ) const;

  //! \f$L^2(\Omega)\f$-inner product with state \c state.
  double_complex
  innerProduct( const Ginla::FVM::State & state ) const;

  //! \f$L^2(\Omega)\f$-norm.
  double
  normalizedScaledL2Norm () const;

  //! Updates the values of the state according to
  //! \f$a \leftarrow \alpha b + \beta a \f$.
  void
  update( const double                  alpha,
          const Ginla::State::Virtual & b,
          const double                  beta
        );

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
  ComplexMultiVector psi_;

  double chi_;

  //! The grid on which the state exists.
  const Teuchos::RCP<const VIO::TpetraMesh::Mesh> mesh_;

};

}

}

#endif // GINLA_FVM_STATE_H
