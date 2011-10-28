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
#ifndef GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_
#define GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <LOCA_Parameter_Vector.H>

// forward declarations
namespace Ginla {
  namespace EpetraFVM {
    class StkMesh;
  }
}

typedef Teuchos::Tuple<double,3> Point;

namespace Ginla {
  namespace MagneticVectorPotential {

class Virtual
{
public:
  Virtual( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh );

  virtual
  ~Virtual();

  //! Sets the parameters in this module.
  //! @return Indicates whether the internal values have changed.
  virtual
  bool
  setParameters( const LOCA::ParameterVector & p ) = 0;

  virtual
  Teuchos::RCP<LOCA::ParameterVector>
  getParameters() const = 0;

  virtual
  Teuchos::RCP<Point>
  getA(const Point & x ) const = 0;

  //! Return the projection of the magnetic vector potential onto the
  //! edge with index \c edgeIndex at the midpoint of the edge.
  virtual
  double
  getAEdgeMidpointProjection( const unsigned int cellIndex,
                              const unsigned int edgeIndex
                            ) const = 0;

protected:
  const Teuchos::RCP<Ginla::EpetraFVM::StkMesh> mesh_;
private:
};

  } // namespace MagneticVectorPotential
} // namespace GL
#endif // GINLA_MAGNETICVECTORPOTENTIAL_VIRTUAL_H_
