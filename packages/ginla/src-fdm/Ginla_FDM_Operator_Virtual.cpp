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

#include "Ginla_FDM_Operator_Virtual.h"

#include "Ginla_Helpers.h"
#include "Ginla_MagneticVectorPotential_Centered.h"
#include "Recti_Grid_Uniform.h"

#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>

// =============================================================================
Ginla::FDM::Operator::Virtual::
Virtual ( const Teuchos::RCP<Recti::Grid::Uniform>                     & grid,
          const Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> & A,
          const Teuchos::RCP<const ComplexMap>                         & domainMap,
          const Teuchos::RCP<const ComplexMap>                         & rangeMap
        ) :
        domainMap_(domainMap),
        rangeMap_(rangeMap),
        grid_ ( grid ),
        A_ ( A ),
        cacheNeedsUpdating_( false ),
        ALeft_( Teuchos::rcp( new RealVector(domainMap) ) ),
        ARight_( Teuchos::rcp( new RealVector(domainMap) ) ),
        AAbove_( Teuchos::rcp( new RealVector(domainMap) ) ),
        ABelow_( Teuchos::rcp( new RealVector(domainMap) ) ),
        dAdH0Left_( Teuchos::rcp( new RealVector(domainMap) ) ),
        dAdH0Right_( Teuchos::rcp( new RealVector(domainMap) ) ),
        dAdH0Above_( Teuchos::rcp( new RealVector(domainMap) ) ),
        dAdH0Below_( Teuchos::rcp( new RealVector(domainMap) ) ),
        firstTime_( true ),
        AB_( Teuchos::null )
{
  // Build the cache for queries to the magnetic vector potential A.
  this->buildACache_();
}
// =============================================================================
Ginla::FDM::Operator::Virtual::~Virtual()
{
}
// =============================================================================
void
Ginla::FDM::Operator::Virtual::setParameters ( const LOCA::ParameterVector & p )
{
    // Set the paramters in A_ and mark the cache obsolete if something changed.
    cacheNeedsUpdating_ = A_->setParameters ( p );
    
    // update the grid
    grid_->updateScaling ( p );
    
    return;
}
// =============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::FDM::Operator::Virtual::getParameters () const
{
  Teuchos::RCP<LOCA::ParameterVector> pA = A_->getParameters();
  Teuchos::RCP<LOCA::ParameterVector> pG = grid_->getParameters();
  
  Teuchos::RCP<LOCA::ParameterVector> pMerge =
      Ginla::Helpers::mergeLocaParameterVectors( *pA, *pG );

  return pMerge;
}
// =============================================================================
Teuchos::RCP<const Recti::Grid::General>
Ginla::FDM::Operator::Virtual::
getGrid() const
{
  return grid_;
}
// =============================================================================
void
Ginla::FDM::Operator::Virtual::
buildACache_() const
{                                     
  Teuchos::ArrayRCP<double> ALeftView  = ALeft_->get1dViewNonConst();
  Teuchos::ArrayRCP<double> ARightView = ARight_->get1dViewNonConst();
  Teuchos::ArrayRCP<double> ABelowView = ABelow_->get1dViewNonConst();
  Teuchos::ArrayRCP<double> AAboveView = AAbove_->get1dViewNonConst();
  
  Teuchos::ArrayRCP<double> dAdH0LeftView  = dAdH0Left_->get1dViewNonConst();
  Teuchos::ArrayRCP<double> dAdH0RightView = dAdH0Right_->get1dViewNonConst();
  Teuchos::ArrayRCP<double> dAdH0BelowView = dAdH0Below_->get1dViewNonConst();
  Teuchos::ArrayRCP<double> dAdH0AboveView = dAdH0Above_->get1dViewNonConst();
  
  // loop over local nodes
  unsigned int localLength = ALeft_->getLocalLength();
  
  for ( unsigned int k=0; k<localLength; k++ )
  {
    int globalIndex = rangeMap_->getGlobalElement ( k );
    
    switch ( grid_->getNodeType ( globalIndex ) )
    {
      case Recti::Grid::Abstract::INTERIOR:
      case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONCAVE:
      case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONCAVE:
      case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONCAVE:
      case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONCAVE:
        ALeftView[k]  = A_->getAx ( *grid_->getXLeft ( globalIndex ) );
        ARightView[k] = A_->getAx ( *grid_->getXRight ( globalIndex ) );
        ABelowView[k] = A_->getAy ( *grid_->getXBelow ( globalIndex ) );
        AAboveView[k] = A_->getAy ( *grid_->getXAbove ( globalIndex ) );
        
        dAdH0LeftView[k]  = A_->getDAxDh0 ( *grid_->getXLeft ( globalIndex ) );
        dAdH0RightView[k] = A_->getDAxDh0 ( *grid_->getXRight ( globalIndex ) );
        dAdH0BelowView[k] = A_->getDAyDh0 ( *grid_->getXBelow ( globalIndex ) );
        dAdH0AboveView[k] = A_->getDAyDh0 ( *grid_->getXAbove ( globalIndex ) );
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONVEX:
        ARightView[k] = A_->getAx ( *grid_->getXRight ( globalIndex ) );
        AAboveView[k] = A_->getAy ( *grid_->getXAbove ( globalIndex ) );
        
        dAdH0RightView[k] = A_->getDAxDh0 ( *grid_->getXRight ( globalIndex ) );
        dAdH0AboveView[k] = A_->getDAyDh0 ( *grid_->getXAbove ( globalIndex ) );
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONVEX:
        ALeftView[k]  = A_->getAx ( *grid_->getXLeft ( globalIndex ) );
        AAboveView[k] = A_->getAy ( *grid_->getXAbove ( globalIndex ) );
        
        dAdH0LeftView[k]  = A_->getDAxDh0 ( *grid_->getXLeft ( globalIndex ) );
        dAdH0AboveView[k] = A_->getDAyDh0 ( *grid_->getXAbove ( globalIndex ) );
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONVEX:
        ALeftView[k]  = A_->getAx ( *grid_->getXLeft ( globalIndex ) );
        ABelowView[k] = A_->getAy ( *grid_->getXBelow ( globalIndex ) );
        
        dAdH0LeftView[k]  = A_->getDAxDh0 ( *grid_->getXLeft ( globalIndex ) );
        dAdH0BelowView[k] = A_->getDAyDh0 ( *grid_->getXBelow ( globalIndex ) );
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONVEX:
        ARightView[k] = A_->getAx ( *grid_->getXRight ( globalIndex ) );
        ABelowView[k] = A_->getAy ( *grid_->getXBelow ( globalIndex ) );
        
        dAdH0RightView[k] = A_->getDAxDh0 ( *grid_->getXRight ( globalIndex ) );
        dAdH0BelowView[k] = A_->getDAyDh0 ( *grid_->getXBelow ( globalIndex ) );
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOM:
        ALeftView[k]  = A_->getAx ( *grid_->getXLeft ( globalIndex ) );
        ARightView[k] = A_->getAx ( *grid_->getXRight ( globalIndex ) );
        AAboveView[k] = A_->getAy ( *grid_->getXAbove ( globalIndex ) );
        
        dAdH0LeftView[k]  = A_->getDAxDh0 ( *grid_->getXLeft ( globalIndex ) );
        dAdH0RightView[k] = A_->getDAxDh0 ( *grid_->getXRight ( globalIndex ) );
        dAdH0AboveView[k] = A_->getDAyDh0 ( *grid_->getXAbove ( globalIndex ) );
        break;

    case Recti::Grid::Abstract::BOUNDARY_RIGHT:
        ALeftView[k]  = A_->getAx ( *grid_->getXLeft ( globalIndex ) );
        ABelowView[k] = A_->getAy ( *grid_->getXBelow ( globalIndex ) );
        AAboveView[k] = A_->getAy ( *grid_->getXAbove ( globalIndex ) );
        
        dAdH0LeftView[k]  = A_->getDAxDh0 ( *grid_->getXLeft ( globalIndex ) );
        dAdH0BelowView[k] = A_->getDAyDh0 ( *grid_->getXBelow ( globalIndex ) );
        dAdH0AboveView[k] = A_->getDAyDh0 ( *grid_->getXAbove ( globalIndex ) );
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOP:
        ALeftView[k]  = A_->getAx ( *grid_->getXLeft ( globalIndex ) );
        ARightView[k] = A_->getAx ( *grid_->getXRight ( globalIndex ) );
        ABelowView[k] = A_->getAy ( *grid_->getXBelow ( globalIndex ) );
        
        dAdH0LeftView[k]  = A_->getDAxDh0 ( *grid_->getXLeft ( globalIndex ) );
        dAdH0RightView[k] = A_->getDAxDh0 ( *grid_->getXRight ( globalIndex ) );
        dAdH0BelowView[k] = A_->getDAyDh0 ( *grid_->getXBelow ( globalIndex ) );
        break;

    case Recti::Grid::Abstract::BOUNDARY_LEFT:
        ARightView[k] = A_->getAx ( *grid_->getXRight ( globalIndex ) );
        ABelowView[k] = A_->getAy ( *grid_->getXBelow ( globalIndex ) );
        AAboveView[k] = A_->getAy ( *grid_->getXAbove ( globalIndex ) );
        
        dAdH0RightView[k] = A_->getDAxDh0 ( *grid_->getXRight ( globalIndex ) );
        dAdH0BelowView[k] = A_->getDAyDh0 ( *grid_->getXBelow ( globalIndex ) );
        dAdH0AboveView[k] = A_->getDAyDh0 ( *grid_->getXAbove ( globalIndex ) );
        break;

    default:
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Illegal not type \"" << grid_->getNodeType ( globalIndex ) << "\"." );
    }
  }
  
  cacheNeedsUpdating_ = false;
}
// =============================================================================