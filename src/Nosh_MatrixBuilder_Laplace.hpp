// @HEADER
//
//    Builder class for the Laplace operator.
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

#ifndef NOSH_MATRIXBUILDER_LAPLACE_H
#define NOSH_MATRIXBUILDER_LAPLACE_H
// =============================================================================
#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  #include <Teuchos_Time.hpp>
#endif
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Epetra_FECrsGraph.h>

#include <stk_mesh/base/Entity.hpp>

#include "Nosh_MatrixBuilder_Virtual.hpp"
// =============================================================================
// forward declarations
namespace Nosh {
class StkMesh;
namespace ScalarField {
class Virtual;
}
}
// =============================================================================
namespace Nosh {
namespace MatrixBuilder {
// =============================================================================
class Laplace: public Virtual
{
public:
Laplace(const Teuchos::RCP<const Nosh::StkMesh> &mesh,
        const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness
        );

// Destructor.
~Laplace();

virtual
const Epetra_Comm &
getComm() const;

virtual
const Epetra_FECrsGraph &
getGraph() const;

virtual
void
apply(const std::map<std::string, double> &params,
      const Epetra_Vector &X,
      Epetra_Vector &Y
      ) const;

virtual
void
applyDKDp(const std::map<std::string, double> &params,
          const std::string & paramName,
          const Epetra_Vector &X,
          Epetra_Vector &Y
          ) const;

virtual
void
fill(Epetra_FECrsMatrix &matrix,
     const std::map<std::string, double> &params
     ) const;

//! Gets the parameter with their initial values.
virtual
const std::map<std::string,double>
getParameters() const;

protected:

private:
const Epetra_FECrsGraph
buildGraph_() const;

void
fill_(Epetra_FECrsMatrix &matrix) const;

void
buildGlobalIndexCache_(const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > &edges) const;

void
buildAlphaCache_(const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > & edges,
                 const Teuchos::ArrayRCP<const double> &edgeCoefficients
                 ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
const Teuchos::RCP<Teuchos::Time> fillTime_;
const Teuchos::RCP<Teuchos::Time> buildLaplaceGraphTime_;
#endif
const Teuchos::RCP<const Nosh::StkMesh> mesh_;
const Teuchos::RCP<const Nosh::ScalarField::Virtual> thickness_;

mutable Teuchos::ArrayRCP<Epetra_IntSerialDenseVector> globalIndexCache_;
mutable bool globalIndexCacheUpToDate_;

const Epetra_FECrsGraph graph_;
mutable Epetra_FECrsMatrix matrixCache_;
mutable bool matrixCacheUpToDate_;

mutable Teuchos::ArrayRCP<double> alphaCache_;
mutable bool alphaCacheUpToDate_;
};
// =============================================================================
} // namespace MatrixBuilder
} // namespace Nosh

#endif // NOSH_MATRIXBUILDER_LAPLACE_H
