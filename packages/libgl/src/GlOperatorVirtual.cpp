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

#include "GlOperatorVirtual.h"

#include <Teuchos_RCP.hpp>

// =============================================================================
GlOperatorVirtual::GlOperatorVirtual ( Teuchos::RCP<GridUniform>             & grid,
                                       Teuchos::RCP<MagneticVectorPotential> & A ) :
        grid_ ( grid ),
        A_ ( A )
{
}
// =============================================================================
GlOperatorVirtual::~GlOperatorVirtual()
{
}
// =============================================================================
void
GlOperatorVirtual::updatePsi ( const Teuchos::RCP<const ComplexVector> psi )
{
    psi_ = psi;
}
// =============================================================================
void
GlOperatorVirtual::setChi ( const double chi )
{
    chi_ = chi;
}
// =============================================================================
void
GlOperatorVirtual::setH0 ( const double h0 )
{
    A_->setH0 ( h0 );
}
// =============================================================================
double
GlOperatorVirtual::getH0 () const
{
    return A_->getH0();
}
// =============================================================================
void
GlOperatorVirtual::setScaling ( const double scaling )
{
    grid_->setScaling( scaling );
    A_->setEdgeLength( scaling );
}
// =============================================================================
double
GlOperatorVirtual::getScaling () const
{
    return grid_->getScaling();
}
// =============================================================================
Teuchos::RCP<const GridVirtual>
GlOperatorVirtual::getGrid() const
{
  return grid_;
}
// =============================================================================