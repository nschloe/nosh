// @HEADER
//
//    Virtual class for matrix constructors.
//    Copyright (C) 2012  Nico Schl\"omer
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

#ifndef NOSH_MATRIXBUILDER_VIRTUAL
#define NOSH_MATRIXBUILDER_VIRTUAL
// =============================================================================
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Epetra_Comm.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>
// =============================================================================
// forward declarations
namespace Nosh {
class StkMesh;
namespace MagneticVectorPotential {
class Virtual;
}
}
// =============================================================================
// =============================================================================
namespace Nosh {
// =============================================================================
namespace MatrixBuilder {
// =============================================================================
class Virtual
{

public:
Virtual();

// Destructor.
~Virtual();

//! Get the underlying communicator.
virtual
const Epetra_Comm &
getComm() const = 0;

//! Get the connectivity graph of the matrix.
virtual
const Epetra_FECrsGraph &
getGraph() const = 0;

//! Y = A(params) * X.
virtual
void
apply(const Teuchos::Array<double> &params,
      const Epetra_Vector &X,
      Epetra_Vector &Y
      ) const = 0;

//! Y = dA/dp(params) * X.
virtual
void
applyDKDp(const Teuchos::Array<double> &params,
          const int paramIndex,
          const Epetra_Vector &X,
          Epetra_Vector &Y
          ) const = 0;

//! Fill a given matrix with the parameter entries as given in params.
virtual
void
fill(Epetra_FECrsMatrix &matrix,
     const Teuchos::Array<double> &params
     ) const = 0;

//! Gets the initial parameters from this module.
virtual
Teuchos::RCP<const Teuchos::Array<double> >
get_p_init() const = 0;

//! Get the parameter names.
virtual
Teuchos::RCP<const Teuchos::Array<std::string> >
get_p_names() const = 0;
};
// =============================================================================
} // namespace MatrixBuilder
} // namespace Nosh

#endif // NOSH_MATRIXBUILDER_VIRTUAL
