/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2009--2010 Nico Schl\"omer, Daniele Avitabile

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

#include "Ginla_Perturbation_Quadrants.h"
#include "Teuchos_Assert.hpp"

// ============================================================================ 
// Constructor
Ginla::Perturbation::Quadrants::
Quadrants( const Teuchos::RCP<const Recti::Grid::Abstract> & grid ):
  Ginla::Perturbation::Virtual(),
  epsilonQuadrant1_(0.0),
  grid_(grid)
{
}
// ============================================================================ 
// Destructor
Ginla::Perturbation::Quadrants::
~Quadrants()
{
}
// ============================================================================ 
// Compute perturbation
std::complex<double>
Ginla::Perturbation::Quadrants::
computePerturbation( int k ) const
{
  Teuchos::RCP<Teuchos::Tuple<double,2> > xy = grid_->getX(k);
  if ( ( (*xy)[0] >= 0.0 ) && ( (*xy)[1] >= 0.0 ) )
      return std::complex<double>(epsilonQuadrant1_,0.0);
  else
      return std::complex<double>(0.0,0.0);
}
// ============================================================================ 
// Set parameters
void Ginla::Perturbation::Quadrants::
setParameters( const LOCA::ParameterVector & p )
{
  TEST_FOR_EXCEPTION ( !p.isParameter ( "Epsilon Quadrant 1" ),
                         std::logic_error,
                         "Label \"Epsilon Quadrant 1\" not valid." );
    epsilonQuadrant1_ = p.getValue ( "Epsilon Quadrant 1" );

  return;
}
// ============================================================================
