/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

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

#ifndef DOMAINCIRCLE_H
#define DOMAINCIRCLE_H

#include "Recti_Domain_Ellipse.h"

namespace Recti
{
  namespace Domain
  {

class Circle:
            public Ellipse
{
public:
    //! Constructor.
    Circle ( double radius );

    //! Destructor.
    virtual
    ~Circle();

protected:
private:
};

  } // namespace Domain
} // namespace Recti

#endif // DOMAINCIRCLE_H
