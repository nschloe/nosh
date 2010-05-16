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

#include "Ginla_Operator_Virtual.h"

#include "Ginla_Helpers.h"
#include "Ginla_MagneticVectorPotential_Centered.h"
#include "Recti_Grid_Uniform.h"

#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>

// =============================================================================
Ginla::Operator::Virtual::
Virtual ( const Teuchos::RCP<Recti::Grid::Uniform>                     & grid,
          const Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> & A,
          const Teuchos::RCP<const ComplexMap>                         & domainMap,
          const Teuchos::RCP<const ComplexMap>                         & rangeMap
          ) :
        domainMap_(domainMap),
        rangeMap_(rangeMap),
        grid_ ( grid ),
        A_ ( A )
{
}
// =============================================================================
Ginla::Operator::Virtual::~Virtual()
{
}
// =============================================================================
void
Ginla::Operator::Virtual::setParameters ( const LOCA::ParameterVector & p )
{
    A_->setParameters ( p );
    grid_->updateScaling ( p );

    return;
}
// =============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::Operator::Virtual::getParameters () const
{
  Teuchos::RCP<LOCA::ParameterVector> pA = A_->getParameters();
  Teuchos::RCP<LOCA::ParameterVector> pG = grid_->getParameters();
  
  Teuchos::RCP<LOCA::ParameterVector> pMerge =
      Ginla::Helpers::mergeLocaParameterVectors( *pA, *pG );

  return pMerge;
}
// =============================================================================
Teuchos::RCP<const Recti::Grid::General>
Ginla::Operator::Virtual::getGrid() const
{
  return grid_;
}
// =============================================================================
