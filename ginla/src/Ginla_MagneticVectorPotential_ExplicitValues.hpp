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
#ifndef GINLA_MAGNETICVECTORPOTENTIAL_EXPLICITVALUES_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_EXPLICITVALUES_H_

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
class ExplicitValues: public Virtual
{
public:
  ExplicitValues( const Teuchos::RCP<Ginla::StkMesh> & mesh,
                  const Teuchos::RCP<const Epetra_MultiVector> & mvp,
                  double mu
                );

  ~ExplicitValues();

  //! Sets the parameters in this module.
  void
  setParameters( const LOCA::ParameterVector & p );

  Teuchos::RCP<LOCA::ParameterVector>
  getParameters() const;

  double
  getAEdgeMidpointProjection( const unsigned int edgeIndex
                            ) const;

  double
  getdAdMuEdgeMidpointProjection( const unsigned int edgeIndex
                                ) const;
  double
  getdAdThetaEdgeMidpointProjection( const unsigned int edgeIndex
                                   ) const;

  double
  getAEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                      const unsigned int edgeIndex
                                    ) const;

  double
  getdAdMuEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                          const unsigned int edgeIndex
                                        ) const;
  double
  getdAdThetaEdgeMidpointProjectionFallback( const unsigned int cellIndex,
                                             const unsigned int edgeIndex
                                           ) const;

protected:
private:
  void
  initializeMvpEdgeMidpointCache_() const;

  void
  initializeMvpEdgeMidpointFallback_() const;

  DoubleVector
  crossProduct_( const DoubleVector u,
                 const DoubleVector v
               ) const;

private:
  const Teuchos::RCP<Ginla::StkMesh> mesh_;
  const Teuchos::RCP<const Epetra_MultiVector> mvp_;
  double mu_;

  Teuchos::ArrayRCP<double> mvpEdgeMidpointProjectionCache_;
  mutable bool mvpEdgeMidpointProjectionCacheUptodate_;
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > mvpEdgeMidpointProjectionCacheFallback_;
  mutable bool mvpEdgeMidpointProjectionCacheFallbackUptodate_;
};
} // namespace MagneticVectorPotential
} // namespace Ginla
#endif // GINLA_MAGNETICVECTORPOTENTIAL_EXPLICITVALUES_H_
