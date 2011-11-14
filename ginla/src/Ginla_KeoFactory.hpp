// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2010, 2011  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER

#ifndef GINLA_KEOFACTORY_H
#define GINLA_KEOFACTORY_H
// =============================================================================
// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include "Ginla_config.h"

#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  #include <Teuchos_Time.hpp>
#endif
#include <Teuchos_Tuple.hpp>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_LinearProblem.h>

#include "Ginla_MagneticVectorPotential.hpp"
// =============================================================================
// forward declarations
namespace Ginla {
    class StkMesh;
}
// =============================================================================
namespace Ginla {
// =============================================================================
class KeoFactory
{
public:
    KeoFactory( const Teuchos::RCP<Ginla::StkMesh>      & mesh,
                const Teuchos::RCP<const Epetra_Vector>            & thickness,
                const Teuchos::RCP<Ginla::MagneticVectorPotential> & mvp
              );

    // Destructor.
    ~KeoFactory();

    void
    updateParameters( const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams
                    ) const;

    const Teuchos::RCP<const LOCA::ParameterVector>
    getMvpParameters() const;

    void
    buildKeo( Epetra_FECrsMatrix & keoMatrix
            ) const;

    const Epetra_FECrsGraph
    buildKeoGraph() const;

protected:
private:
#ifdef GINLA_TEUCHOS_TIME_MONITOR
    const Teuchos::RCP<Teuchos::Time> buildKeoTime_;
#endif
    const Teuchos::RCP<Ginla::StkMesh> mesh_;
    const Teuchos::RCP<const Epetra_Vector>       thickness_;
    const Teuchos::RCP<Ginla::MagneticVectorPotential> mvp_;
};
// =============================================================================
} // namespace Ginla

#endif // GINLA_KEOFACTORY_H
