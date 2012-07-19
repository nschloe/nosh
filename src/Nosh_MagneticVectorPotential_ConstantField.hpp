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
#ifndef NOSH_MAGNETICVECTORPOTENTIAL_CONSTANTFIELD_H_
#define NOSH_MAGNETICVECTORPOTENTIAL_CONSTANTFIELD_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Epetra_MultiVector.h>
#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

#include "Nosh_MagneticVectorPotential_Virtual.hpp"
#include "Nosh_StkMesh.hpp"

typedef Teuchos::SerialDenseVector<int,double> DoubleVector;

namespace Nosh {
namespace MagneticVectorPotential {
class ConstantField : public Virtual
{
public:
ConstantField(const Teuchos::RCP<Nosh::StkMesh> &mesh,
              const Teuchos::RCP<DoubleVector> &b,
              const Teuchos::RCP<DoubleVector> &u = Teuchos::null
              );

virtual
~ConstantField();

//! Get the parameter names.
virtual
Teuchos::RCP<const Teuchos::Array<std::string> >
get_p_names() const;

//! Gets the current parameters from this module.
virtual
Teuchos::RCP<const Teuchos::Array<double> >
get_p_init() const;

virtual
double
getAEdgeMidpointProjection(const unsigned int edgeIndex,
                           const Teuchos::Array<double> &mvpParams
                           ) const;

virtual
double
getdAdPEdgeMidpointProjection(const unsigned int edgeIndex,
                              const Teuchos::Array<double> &mvpParams,
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

void
rotate_(DoubleVector &v,
        const DoubleVector &u,
        const double theta
        ) const;

void
dRotateDTheta_(DoubleVector &v,
               const DoubleVector &u,
               const double theta
               ) const;

DoubleVector
crossProduct_(const DoubleVector u,
              const DoubleVector v
              ) const;

private:
const Teuchos::RCP<Nosh::StkMesh> mesh_;
const Teuchos::RCP<const DoubleVector> b_;
const Teuchos::RCP<const DoubleVector> u_;
mutable DoubleVector rotatedBCache_;
mutable double rotatedBCacheAngle_;
mutable DoubleVector dRotatedBDThetaCache_;
mutable double rotateddBdThetaCacheAngle_;

Teuchos::ArrayRCP<DoubleVector> edgeCache_;
mutable bool edgeCacheUptodate_;
};
} // namespace MagneticVectorPotential
} // namespace Nosh
#endif // NOSH_MAGNETICVECTORPOTENTIAL_CONSTANTFIELD_H_
