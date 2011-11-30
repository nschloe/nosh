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
#include <Epetra_MultiVector.h>
#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

#include "Ginla_StkMesh.hpp"

namespace Ginla {

class MagneticVectorPotential
{
public:
  MagneticVectorPotential( const Teuchos::RCP<Ginla::StkMesh> & mesh,
                           const Teuchos::RCP<const Epetra_MultiVector>  & mvp,
                           double mu
                         );

  ~MagneticVectorPotential();

  //! Sets the parameters in this module.
  void
  setParameters( const LOCA::ParameterVector & p );

  Teuchos::RCP<LOCA::ParameterVector>
  getParameters() const;

  Teuchos::RCP<Point>
  getA(const Point & x ) const;

  double
  getAEdgeMidpointProjection( const unsigned int cellIndex,
                              const unsigned int edgeIndex
                            ) const;

  double
  getdAdMuEdgeMidpointProjection( const unsigned int cellIndex,
                                  const unsigned int edgeIndex
                                ) const;

protected:
private:
  void
  initializeEdgeMidpointProjectionCache_() const;

private:
  const Teuchos::RCP<Ginla::StkMesh> mesh_;
  const Teuchos::RCP<const Epetra_MultiVector> mvp_;
  double mu_;

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeMidpointProjectionCache_;
  mutable bool edgeMidpointProjectionCacheUpToDate_;

};
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_H_
