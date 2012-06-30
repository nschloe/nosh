// @HEADER
//
//    Query routines for the magnetic vector potential.
//    Copyright (C) 2011, 2012  Nico Schl\"omer
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
#ifndef GINLA_MAGNETICVECTORPOTENTIAL_CONSTANTINSPACE_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_CONSTANTINSPACE_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Epetra_MultiVector.h>
#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

#include "Ginla_MagneticVectorPotential_Virtual.hpp"
#include "Ginla_StkMesh.hpp"

typedef Teuchos::SerialDenseVector<int,double> DoubleVector;

namespace Ginla {
namespace MagneticVectorPotential {
class ConstantInSpace : public Virtual
{
public:
ConstantInSpace( const Teuchos::RCP<Ginla::StkMesh> &mesh,
                 const Teuchos::RCP<DoubleVector> &p,
                 double mu,
                 double theta = 0.0,
                 const Teuchos::RCP<DoubleVector> &u = Teuchos::null
                 );

~ConstantInSpace();

//! Sets the parameters in this module.
void
setParameters( const LOCA::ParameterVector &p );

Teuchos::RCP<LOCA::ParameterVector>
getParameters() const;

double
getAEdgeMidpointProjection(const unsigned int edgeIndex
                           ) const;

double
getdAdPEdgeMidpointProjection(const unsigned int edgeIndex,
                              const unsigned int parameterIndex
                              ) const;

protected:
private:
DoubleVector
getRawA_( const DoubleVector &x ) const;

DoubleVector
getRawDADTheta_( const DoubleVector &x ) const;

void
initializeEdgeCache_() const;

DoubleVector
rotate_( const DoubleVector &v,
         const DoubleVector &u,
         const double sinTheta,
         const double cosTheta
         ) const;

DoubleVector
dRotateDTheta_( const DoubleVector &v,
                const DoubleVector &u,
                const double sinTheta,
                const double cosTheta
                ) const;

DoubleVector
crossProduct_( const DoubleVector u,
               const DoubleVector v
               ) const;

private:
const Teuchos::RCP<Ginla::StkMesh> mesh_;
const Teuchos::RCP<const DoubleVector> b_;
DoubleVector rotatedB_;
DoubleVector dRotatedBDTheta_;
double mu_;
double theta_;
const Teuchos::RCP<const DoubleVector> u_;

Teuchos::ArrayRCP<DoubleVector> edgeCache_;
mutable bool edgeCacheUptodate_;
};
} // namespace MagneticVectorPotential
} // namespace Ginla
#endif // GINLA_MAGNETICVECTORPOTENTIAL_CONSTANTINSPACE_H_
