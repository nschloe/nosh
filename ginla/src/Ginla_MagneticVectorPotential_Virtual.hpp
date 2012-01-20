// @HEADER
//
//    Query routines for the magnetic vector potential.
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
#ifndef GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_
// =============================================================================
#include <LOCA_Parameter_Vector.H>
#include <Teuchos_RCP.hpp>
// =============================================================================
namespace Ginla {
namespace MagneticVectorPotential {
class Virtual
{
public:
Virtual();

~Virtual();

//! Sets the parameters in this module.
virtual void
setParameters( const LOCA::ParameterVector &p ) = 0;

virtual Teuchos::RCP<LOCA::ParameterVector>
getParameters() const = 0;

virtual double
getAEdgeMidpointProjection( const unsigned int edgeIndex
                            ) const = 0;

virtual double
getdAdMuEdgeMidpointProjection( const unsigned int edgeIndex
                                ) const = 0;
virtual double
getdAdThetaEdgeMidpointProjection( const unsigned int edgeIndex
                                   ) const = 0;

virtual double
getAEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                    const unsigned int edgeIndex
                                    ) const = 0;

virtual double
getdAdMuEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                        const unsigned int edgeIndex
                                        ) const = 0;
virtual double
getdAdThetaEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                           const unsigned int edgeIndex
                                           ) const = 0;

protected:
private:
};
} // namespace MagneticVectorPotential
} // namespace Ginla
#endif // GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_
