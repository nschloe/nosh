/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2011  Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef GINLA_MAGNETICVECTORPOTENTIAL_CUSTOM_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_CUSTOM_H_

#include "Ginla_MagneticVectorPotential_Virtual.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Epetra_MultiVector.h>
#include <Teuchos_Array.hpp>
#include <LOCA_Parameter_Vector.H>

namespace Ginla {
  namespace MagneticVectorPotential {

class Custom:
  public Virtual
{
public:
  Custom( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
          const Teuchos::RCP<const Epetra_MultiVector>  & mvp,
          double mu
        );

  virtual
  ~Custom();

  //! Sets the parameters in this module.
  //! @return Indicates whether the internal values have changed.
  virtual
  bool
  setParameters( const LOCA::ParameterVector & p );

  virtual
  Teuchos::RCP<LOCA::ParameterVector>
  getParameters() const;

  virtual
  Teuchos::RCP<Point>
  getA(const Point & x ) const;

  virtual
  double
  getAEdgeMidpointProjection( const unsigned int cellIndex,
                              const unsigned int edgeIndex
                            ) const;

protected:
private:

  void
  initializeEdgeMidpointProjectionCache_() const;

private:

  double mu_;

  const Teuchos::RCP<const Epetra_MultiVector> mvp_;

  const Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > edgeMidpointProjectionCache_;
  mutable bool edgeMidpointProjectionCacheUpToDate_;

};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_CUSTOM_H_
