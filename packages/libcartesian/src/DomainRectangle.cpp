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

#include "DomainRectangle.h"

// ============================================================================
DomainRectangle::DomainRectangle ( const double a,
                                   const double b ) :
        DomainVirtual(),
        a_(a),
        b_(b)
{
}
// ============================================================================
DomainRectangle::~DomainRectangle()
{
}
// ============================================================================
Teuchos::Tuple<double,4>
DomainRectangle::getBoundingBox () const
{
    return Teuchos::tuple( 0.0, 0.0, a_, b_ );
}
// ============================================================================
bool
DomainRectangle::isInDomain ( const DoubleTuple & x ) const
{
    return    x[0]>=-tolerance_ && x[0]<=a_+tolerance_
           && x[1]>=-tolerance_ && x[1]<=b_+tolerance_;
}
// ============================================================================
