// @HEADER
//
//    Container class that hosts the kinetic energy operator.
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

#ifndef GINLA_KEOCONTAINER_H
#define GINLA_KEOCONTAINER_H
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
#include <Teuchos_Array.hpp>
#include <Teuchos_SerialDenseVector.hpp>

#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <stk_mesh/base/Entity.hpp>
// =============================================================================
typedef Teuchos::SerialDenseVector<int,double> DoubleVector;
// =============================================================================
// forward declarations
namespace Ginla {
class StkMesh;
namespace MagneticVectorPotential {
class Virtual;
}
}
// =============================================================================
// =============================================================================
namespace Ginla {
// =============================================================================
class KeoContainer
{
public:
enum EMatrixType { MATRIX_TYPE_REGULAR,
                   MATRIX_TYPE_DMU,
                   MATRIX_TYPE_DTHETA};

public:
KeoContainer(const Teuchos::RCP<const Ginla::StkMesh> &mesh,
             const Teuchos::RCP<const Epetra_Vector> &thickness,
             const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> &mvp
             );

// Destructor.
~KeoContainer();

const Epetra_Comm &
getComm() const;

//void
//updateParameters(const Teuchos::Array<double> &mvpParams
//                 ) const;

Teuchos::RCP<const Epetra_FECrsMatrix>
getKeo(const Teuchos::Array<double> & mvpParams) const;

Teuchos::RCP<const Epetra_FECrsMatrix>
getKeoDp(const int paramIndex,
         const Teuchos::Array<double> & mvpParams
         ) const;

protected:

private:
const Teuchos::RCP<Epetra_FECrsGraph>
buildKeoGraph_() const;

void
fillKeo_( const Teuchos::RCP<Epetra_FECrsMatrix> &keoMatrix,
          const Teuchos::Array<double> &mvpParams,
          void (KeoContainer::*filler)(const int, const Teuchos::Array<double>&, double*) const
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
#ifdef GINLA_TEUCHOS_TIME_MONITOR
const Teuchos::RCP<Teuchos::Time> keoFillTime_;
const Teuchos::RCP<Teuchos::Time> buildKeoGraphTime_;
#endif
const Teuchos::RCP<const Ginla::StkMesh> mesh_;
const Teuchos::RCP<const Epetra_Vector> thickness_;
const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> mvp_;

mutable Teuchos::ArrayRCP<Epetra_IntSerialDenseVector> globalIndexCache_;
mutable bool globalIndexCacheUpToDate_;

const Teuchos::RCP<const Epetra_FECrsGraph> keoGraph_;
const Teuchos::RCP<Epetra_FECrsMatrix> keo_;
mutable Teuchos::Array<double> keoBuildParameters_;
mutable Teuchos::Array<Teuchos::RCP<Epetra_FECrsMatrix> > keoDp_;
mutable Teuchos::Array<Teuchos::Array<double> > keoDpBuildParameters_;

mutable Teuchos::ArrayRCP<double> alphaCache_;
mutable bool alphaCacheUpToDate_;
mutable unsigned int paramIndex_;
};
// =============================================================================
} // namespace Ginla

#endif // GINLA_KEOCONTAINER_H
