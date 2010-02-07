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

#ifndef GRIDUNIFORM_H
#define GRIDUNIFORM_H

#include "Grid.h"

class GridUniform:
            virtual public Grid
{
public:
    GridUniform ( const Teuchos::RCP<const DomainVirtual> & domain,
                  const double                              h,
                  const double                              scaling
                );

    double
    getUniformH() const;

private:
    double h_;

};

#endif // GRIDUNIFORM_H
