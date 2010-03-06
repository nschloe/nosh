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

#include "GL_Operator_Virtual.h"

#include "GL_Helpers.h"

#include <Teuchos_RCP.hpp>

// =============================================================================
GL::Operator::Virtual::Virtual ( Teuchos::RCP<GridUniform>             & grid,
                                 Teuchos::RCP<GL::MagneticVectorPotential::Centered> & A ) :
        psi_ ( Teuchos::null ),
        chi_ ( 0.0 ),
        grid_ ( grid ),
        A_ ( A )
{
}
// =============================================================================
GL::Operator::Virtual::~Virtual()
{
}
// =============================================================================
void
GL::Operator::Virtual::updatePsi ( const Teuchos::RCP<const ComplexVector> psi )
{
    psi_ = psi;
}
// =============================================================================
void
GL::Operator::Virtual::setParameters ( const LOCA::ParameterVector & p )
{
    chi_ = p.getValue ( "chi" );

    A_->setParameters ( p );
    grid_->updateScaling ( p );
    
    return;
}
// =============================================================================
Teuchos::RCP<LOCA::ParameterVector>
GL::Operator::Virtual::getParameters () const
{
  Teuchos::RCP<LOCA::ParameterVector> pA = A_->getParameters();
  Teuchos::RCP<LOCA::ParameterVector> pG = grid_->getParameters();
  
  Teuchos::RCP<LOCA::ParameterVector> pMerge =
      GL::Helpers::mergeLocaParameterVectors( *pA, *pG );
      
  pMerge->addParameter( "chi", chi_ );

  return pMerge;
}
// =============================================================================
Teuchos::RCP<const Grid>
GL::Operator::Virtual::getGrid() const
{
  return grid_;
}
// =============================================================================
