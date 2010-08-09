/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl"omer

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

#ifndef GINLA_STATE_VIRTUAL_H
#define GINLA_STATE_VIRTUAL_H

#include "Ginla_Typedefs.h"

namespace LOCA {
  class ParameterVector;
}

namespace Ginla {
namespace State {
  class Virtual
{
public:

  //! Constructor.
  Virtual();

  virtual
  ~Virtual();
  
  //! Save the state to file \c fileName together with the parameters \c p.
  virtual void
  save( const std::string            & fileName,
        const Teuchos::ParameterList & p
      ) const = 0;
  
  //! Just plain save the file to \c fileName.
  virtual void
  save( const std::string & fileName
      ) const = 0;
      
  virtual
  double
  freeEnergy() const = 0;
  
  virtual
  double
  normalizedScaledL2Norm() const = 0;
  
  virtual
  int
  getVorticity() const = 0;
  
  virtual
  void
  update( const double                  alpha,
          const Ginla::State::Virtual & b,
          const double                  beta
        ) = 0;

  virtual
  Teuchos::RCP<const ComplexVector>
  getPsi () const = 0;

  virtual
  Teuchos::RCP<ComplexVector>
  getPsiNonConst () = 0;
  
  virtual
  double
  getChi () const = 0;


protected:
private:
  
};

}
}

#endif // GINLA_STATE_VIRTUAL_H
