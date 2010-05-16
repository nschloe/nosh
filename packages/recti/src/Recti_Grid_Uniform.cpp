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

#include "Recti_Grid_Uniform.h"

// =============================================================================
Recti::Grid::Uniform::
Uniform ( const Teuchos::RCP<const Domain::Abstract> & domain,
          const double                                  h
        ) :
        General ( domain, Teuchos::tuple ( h,h ) )
{
}
// =============================================================================
Recti::Grid::Uniform::
Uniform ( const double                h,
          const UIntTuple           & numCells,
          const Teuchos::Array<int> & kBB,
          const Teuchos::Array<int> & boundaryNodes,
          const double                scaling,
          const DoubleTuple         & origin
        ) :
        General ( Teuchos::tuple ( h,h ), numCells, kBB, boundaryNodes, scaling, origin )
{
}
// =============================================================================
double
Recti::Grid::Uniform::
getUniformH() const
{
    return h_[0];
}
// ============================================================================