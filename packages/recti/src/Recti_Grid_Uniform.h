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

#ifndef GRIDUNIFORM_H
#define GRIDUNIFORM_H

#include "Recti_Grid_General.h"

namespace Recti
{
  namespace Grid
  {

class Uniform:
            virtual public General
{
public:
    Uniform ( const Teuchos::RCP<const Domain::Abstract> & domain,
              const double                                 h
            );

    Uniform ( const double                h,
              const UIntTuple           & numCells,
              const Teuchos::Array<int> & kBB,
              const Teuchos::Array<int> & boundaryNodes,
              const double                scaling,
              const DoubleTuple         & origin
            );

    double
    getUniformH() const;

private:
};
  } // namespace Grid
} // namespace Recti

#endif // GRIDUNIFORM_H
