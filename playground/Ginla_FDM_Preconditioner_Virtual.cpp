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

#include "Ginla_Preconditioner_Virtual.h"

// =============================================================================
Ginla::Preconditioner::Virtual::
Virtual ( const Teuchos::RCP<Recti::Grid::Uniform> & grid,
          const Teuchos::RCP<const ComplexMap>     & domainMap,
          const Teuchos::RCP<const ComplexMap>     & rangeMap
        ) :
        domainMap_(domainMap),
        rangeMap_(rangeMap),
        grid_ ( grid ),
        firstTime_( true ),
        AB_( Teuchos::null )
{
}
// =============================================================================
Ginla::Preconditioner::Virtual::
~Virtual()
{
}
// =============================================================================


