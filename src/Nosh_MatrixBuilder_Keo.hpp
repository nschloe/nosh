// @HEADER
//
//    Builder class for the kinetic energy operator.
//    Copyright (C) 2010--2012  Nico Schl\"omer
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

#ifndef NOSH_MATRIXBUILDER_KEO_H
#define NOSH_MATRIXBUILDER_KEO_H
// =============================================================================
#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  #include <Teuchos_Time.hpp>
#endif
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Epetra_FECrsGraph.h>

#include <stk_mesh/base/Entity.hpp>

#include "Nosh_MatrixBuilder_Virtual.hpp"
// =============================================================================
typedef Teuchos::SerialDenseVector<int,double> DoubleVector;
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
namespace MatrixBuilder {
// =============================================================================
class Keo: public Virtual
{
public:
enum EMatrixType { MATRIX_TYPE_REGULAR,
                   MATRIX_TYPE_DMU,
                   MATRIX_TYPE_DTHETA};

public:
Keo(const Teuchos::RCP<const Nosh::StkMesh> &mesh,
    const Teuchos::RCP<const Epetra_Vector> &thickness,
    const Teuchos::RCP<const Nosh::MagneticVectorPotential::Virtual> &mvp
    );

// Destructor.
~Keo();

const Epetra_Comm &
getComm() const;

const Epetra_FECrsGraph &
getGraph() const;

void
apply(const Teuchos::Array<double> &mvpParams,
      const Epetra_Vector &X,
         Epetra_Vector &Y
         ) const;

void
applyDKDp(const Teuchos::Array<double> &mvpParams,
          const int paramIndex,
          const Epetra_Vector &X,
          Epetra_Vector &Y
          ) const;

void
fill(Epetra_FECrsMatrix &matrix,
     const Teuchos::Array<double> &mvpParams
     ) const;

protected:

private:
const Epetra_FECrsGraph
buildKeoGraph_() const;

void
fillKeo_( Epetra_FECrsMatrix &keoMatrix,
          const Teuchos::Array<double> &mvpParams,
          void (Keo::*filler)(const int, const Teuchos::Array<double>&, double*) const
          ) const;

void
fillerRegular_(const int k,
               const Teuchos::Array<double> &mvpParams,
               double * v
               ) const;

void
fillerDp_(const int k,
          const Teuchos::Array<double> &mvpParams,
          double * v
          ) const;

void
buildGlobalIndexCache_( const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > &edges ) const;

void
buildAlphaCache_( const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> > & edges,
                  const Teuchos::ArrayRCP<const double> &edgeCoefficients
                ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
const Teuchos::RCP<Teuchos::Time> keoFillTime_;
const Teuchos::RCP<Teuchos::Time> buildKeoGraphTime_;
#endif
const Teuchos::RCP<const Nosh::StkMesh> mesh_;
const Teuchos::RCP<const Epetra_Vector> thickness_;
const Teuchos::RCP<const Nosh::MagneticVectorPotential::Virtual> mvp_;

mutable Teuchos::ArrayRCP<Epetra_IntSerialDenseVector> globalIndexCache_;
mutable bool globalIndexCacheUpToDate_;

const Epetra_FECrsGraph keoGraph_;
mutable Epetra_FECrsMatrix keoCache_;
mutable Teuchos::Array<double> keoBuildParameters_;
mutable Epetra_FECrsMatrix keoDpCache_;


mutable Teuchos::ArrayRCP<double> alphaCache_;
mutable bool alphaCacheUpToDate_;
mutable unsigned int paramIndex_;
};
// =============================================================================
} // namespace MatrixBuilder
} // namespace Nosh

#endif // NOSH_MATRIXBUILDER_KEO_H
