/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

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

#ifndef GINLA_PARAMETERHOST_VIRTUAL_H
#define GINLA_PARAMETERHOST_VIRTUAL_H

#include <Teuchos_RCP.hpp>

namespace LOCA {
  class ParameterVector;
}

namespace Ginla {

namespace ParameterHost {

class Virtual
{
  public:
  Virtual();
  
  virtual
  ~Virtual();
  
  virtual
  Teuchos::RCP<LOCA::ParameterVector>
  getParameters() const = 0;
  
  protected:
  private:
};

}

}

#endif // GINLA_PARAMETERHOST_VIRTUAL_H
