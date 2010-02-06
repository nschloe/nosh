/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl\"omer

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

#include "DomainPolygon.h"

// ============================================================================
DomainPolygon::DomainPolygon ( const Teuchos::Array<Teuchos::Tuple<double,2> > & polygonPoints,
                               const Teuchos::Tuple<double,2>                  & anchorPoint
                             ) :
        DomainVirtual (),
        polygonPoints_ ( polygonPoints )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// ============================================================================
DomainPolygon::~DomainPolygon()
{
}
// ============================================================================
bool
DomainPolygon::isInDomain ( const Teuchos::Array<double> & x ) const
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
    return true;
}
// ============================================================================
