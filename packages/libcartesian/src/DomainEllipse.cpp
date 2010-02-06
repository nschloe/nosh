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

#include "DomainEllipse.h"

// ============================================================================
DomainEllipse::DomainEllipse ( double a, double b ) :
        a_(a),
        b_(b)
{
}
// ============================================================================
DomainEllipse::~DomainEllipse()
{
}
// ============================================================================
Teuchos::Tuple<double,4>
DomainEllipse::getBoundingBox () const
{
    return Teuchos::tuple( -a_, -b_, a_, b_ );
}
// ============================================================================
bool
DomainEllipse::isInDomain ( const DoubleTuple & x ) const
{
    return    x[0]*x[0]/(a_*a_) + x[1]*x[1] / (b_*b_) <= 1.0;
}
// ============================================================================