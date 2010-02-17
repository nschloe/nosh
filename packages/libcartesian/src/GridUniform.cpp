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

#include "GridUniform.h"

// =============================================================================
GridUniform::GridUniform ( const Teuchos::RCP<const DomainVirtual> & domain,
                           const double                              h
                         ) :
        Grid ( domain, Teuchos::tuple ( h,h ) )
{
}
// =============================================================================
GridUniform::GridUniform ( const double                h,
                           const UIntTuple           & numCells,
                           const Teuchos::Array<int> & kBB,
                           const Teuchos::Array<int> & boundaryNodes,
                           const double                scaling,
                           const DoubleTuple         & origin
                         ) :
        Grid ( Teuchos::tuple ( h,h ), numCells, kBB, boundaryNodes, scaling, origin )
{
}
// =============================================================================
double
GridUniform::getUniformH() const
{
    return h_[0];
}
// ============================================================================