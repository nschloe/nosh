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

#ifndef GINLA_EPETRAFVM_STATE_H
#define GINLA_EPETRAFVM_STATE_H
// =============================================================================
// includes
// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
// =============================================================================
// forward declarations
namespace LOCA {
  class ParameterVector;
}
namespace Ginla {
  namespace EpetraFVM {
    class StkMesh;
  }
}
namespace stk {
  namespace mesh {
    class MetaData;
    class BulkData;
  }
}
// =============================================================================
namespace Ginla {
namespace EpetraFVM {
    class State
{
public:

  //! Constructor.
  State( const Epetra_Vector                                 & psi,
         const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh
       );

  //! Constructor without \f$\psi\f$. The values will be initialized to 0.
  State( const Teuchos::RCP<const Epetra_Map>                & map,
         const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh
       );

  //! Constructor solely with comminicator and grid. The values will be initialized to 0.
  State( const Teuchos::RCP<const Epetra_Comm>               & comm,
         const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh
       );

  //! Const getter.
  Teuchos::RCP<const Epetra_Vector>
  getPsi () const;

  //! Nonconst getter for the values.
  Teuchos::RCP<Epetra_Vector>
  getPsiNonConst ();

  const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh>
  getMesh () const;

  //! Save the state to file \c fileName together with the parameters \c p.
  void
  save( const int                      index,
        const Teuchos::ParameterList & p
      ) const;

  //! Just plain save the file to \c fileName.
  void
  save( const int index
      ) const;

  //! \f$L^2(\Omega)\f$-inner product with state \c state.
  double
  innerProduct( const Ginla::EpetraFVM::State & state ) const;

  //! \f$L^2(\Omega)\f$-norm.
  double
  normalizedScaledL2Norm () const;

  //! Updates the values of the state according to
  //! \f$a \leftarrow \alpha b + \beta a \f$.
  void
  update( const double                    alpha,
          const Ginla::EpetraFVM::State & b,
          const double                    beta
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
  Epetra_Vector psi_;

  //! The grid on which the state exists.
  const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> mesh_;

private:

  void
  mergePsi_( const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
             const Epetra_Vector                                 & psi
           ) const;

};

}

}

#endif // GINLA_EPETRAFVM_STATE_H
