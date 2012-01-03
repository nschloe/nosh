// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2011  Nico Schl\"omer
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
#ifndef GINLA_MAGNETICVECTORPOTENTIAL_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Epetra_MultiVector.h>
#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

#include "Ginla_StkMesh.hpp"

typedef Teuchos::SerialDenseVector<int,double> DoubleVector;

namespace Ginla {

class MagneticVectorPotential
{
public:
  MagneticVectorPotential( const Teuchos::RCP<Ginla::StkMesh> & mesh,
                           const Teuchos::RCP<const Epetra_MultiVector>  & mvp,
                           double mu = 0.0,
                           double theta = 0.0,
                           const DoubleVector u = DoubleVector(3)
                         );

  ~MagneticVectorPotential();

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
  rotate_( const DoubleVector & v,
           const DoubleVector & u,
           const double sinTheta,
           const double cosTheta
         ) const;

  DoubleVector
  crossProduct_( const DoubleVector u,
                 const DoubleVector v
               ) const;

private:
  const Teuchos::RCP<Ginla::StkMesh> mesh_;
  const Teuchos::RCP<const Epetra_MultiVector> mvp_;
  double mu_;
  double theta_;
  double sinTheta_;
  double cosTheta_;
  DoubleVector u_;

  Teuchos::ArrayRCP<DoubleVector> mvpEdgeMidpoint_;
  Teuchos::ArrayRCP<DoubleVector> edges_;
  mutable bool mvpEdgeMidpointUpToDate_;

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> > mvpEdgeMidpointFallback_;
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<DoubleVector> > edgesFallback_;
  mutable bool mvpEdgeMidpointFallbackUpToDate_;

};
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_H_
