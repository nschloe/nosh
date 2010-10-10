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

#include "Ginla_FDM_Operator_BCCentral.h"

#include <Teuchos_Array.hpp>

#include "Recti_Grid_Uniform.h"
#include "Ginla_Komplex_DoubleMatrix.h"

// =============================================================================
Ginla::FDM::Operator::BCCentral::
BCCentral ( const Teuchos::RCP<Recti::Grid::Uniform>                    & grid,
            const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & A,
            const Teuchos::RCP<const ComplexMap>                        & domainMap,
            const Teuchos::RCP<const ComplexMap>                        & rangeMap
          ) :
        Ginla::FDM::Operator::Virtual ( grid, A, domainMap, rangeMap )
{
}
// =============================================================================
Ginla::FDM::Operator::BCCentral::
~BCCentral()
{
}
// =============================================================================
Teuchos::RCP<Ginla::State::Virtual>
Ginla::FDM::Operator::BCCentral::
getF( const Teuchos::RCP<const Ginla::FDM::State> & state ) const
{ 
  // initialize F
  Teuchos::RCP<Ginla::State::Virtual> F =
      Teuchos::rcp( new Ginla::FDM::State( state->getPsi()->getMap(),
                                           state->getGrid() ) );
                                      
  Teuchos::ArrayRCP<double_complex> FView = F->getPsiNonConst()->get1dViewNonConst();
  
  // rebuild cache for A?
  if ( cacheNeedsUpdating_ )
    this->buildACache_();
  
  Teuchos::ArrayRCP<const double> ALeftView  = ALeft_->get1dView();
  Teuchos::ArrayRCP<const double> ARightView = ARight_->get1dView();
  Teuchos::ArrayRCP<const double> ABelowView = ABelow_->get1dView();
  Teuchos::ArrayRCP<const double> AAboveView = AAbove_->get1dView();
  
  Teuchos::ArrayRCP<const double_complex> psiView = state->getPsi()->get1dView();
  double chi = state->getChi();
  
  double h = grid_->getUniformH();
  
  // loop over the nodes
  unsigned int localLength = F->getPsi()->getLocalLength();
  
  for ( unsigned int localIndex=0; localIndex<localLength; localIndex++ )
  {
      int globalIndex = rangeMap_->getGlobalElement ( localIndex );
      FView[localIndex] = this->getFEntry_ ( psiView,
                                             chi,
                                             localIndex,
                                             globalIndex,
                                             ALeftView,
                                             ARightView,
                                             ABelowView,
                                             AAboveView,
                                             h );
  }
  
  return F;
}
// =============================================================================
double_complex
Ginla::FDM::Operator::BCCentral::
getFEntry_ ( Teuchos::ArrayRCP<const double_complex> & psiView,
             const double                              chi,
             const int localIndex,
             const int globalIndex,
             Teuchos::ArrayRCP<const double> & ALeftView,
             Teuchos::ArrayRCP<const double> & ARightView,
             Teuchos::ArrayRCP<const double> & ABelowView,
             Teuchos::ArrayRCP<const double> & AAboveView,
             const double                      h
           ) const
{
    double_complex res;
    double_complex psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;
      
    Recti::Grid::Abstract::nodeType nt = grid_->getNodeType ( globalIndex );
    
    psiK = psiView[ localIndex ];
    
    switch ( nt )
    {
    case Recti::Grid::Abstract::INTERIOR:
        // TODO Gets the local index. ==> Only works on one core.
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( psiK* ( -4.0 )
                + psiKLeft*  exp ( IM*ALeftView[localIndex] *h ) + psiKRight* exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow* exp ( IM*ABelowView[localIndex]*h ) + psiKAbove* exp ( -IM*AAboveView[localIndex]*h ) ) / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
              
        res *= exp ( IM*chi );
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONVEX:
        // -------------------------------------------------------------------
        // interior equation, then outward derivative substituted
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiK      * 4.0
                + psiKRight * 2.0 * exp ( -IM*ARightView[localIndex]*h )
                + psiKAbove * 2.0 * exp ( -IM*AAboveView[localIndex]*h ) ) / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONVEX:
        // -----------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiK      * 4.0
                + psiKLeft  * 2.0 * exp ( IM*ALeftView[localIndex] *h )
                + psiKAbove * 2.0 * exp ( -IM*AAboveView[localIndex]*h ) ) / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
        res *= exp ( IM*chi );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONVEX:
        // -----------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];

        res = ( - psiK      * 4.0
                + psiKLeft  * 2.0 * exp ( IM*ALeftView[localIndex] *h )
                + psiKBelow * 2.0 * exp ( IM*ABelowView[localIndex]*h ) ) / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
        res *= exp ( IM*chi );
        // -----------------------------------------------------------------------

        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONVEX:
        // -----------------------------------------------------------------------
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];

        res = ( - psiK      * 4.0
                + psiKRight * 2.0 * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow * 2.0 * exp ( IM*ABelowView[localIndex]*h ) ) / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
        res *= exp ( IM*chi );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOM:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiK      * 4.0
                + psiKLeft  *       exp ( IM*ALeftView[localIndex] *h ) + psiKRight       * exp ( -IM*ARightView[localIndex]*h )
                + psiKAbove * 2.0 * exp ( -IM*AAboveView[localIndex]*h ) )
              / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_RIGHT:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiK      * 4.0
                + psiKLeft  * 2.0 * exp ( IM*ALeftView[localIndex] *h )
                + psiKBelow       * exp ( IM*ABelowView[localIndex]*h ) + psiKAbove * exp ( -IM*AAboveView[localIndex]*h ) )
              / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOP:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];

        res = ( - psiK      * 4.0
                + psiKLeft        * exp ( IM*ALeftView[localIndex] *h ) + psiKRight * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow * 2.0 * exp ( IM*ABelowView[localIndex]*h ) )
              / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_LEFT:
        // -------------------------------------------------------------------
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiK      * 4.0
                + psiKRight * 2.0 * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow * exp ( IM*AAboveView[localIndex]*h ) + psiKAbove       * exp ( -IM*AAboveView[localIndex]*h ) )
              / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONCAVE:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( + psiKLeft  * exp ( IM*ALeftView[localIndex] *h )
                - psiKRight * exp ( -IM*ARightView[localIndex]*h )
                - psiKBelow * exp ( IM*ABelowView[localIndex]*h )
                + psiKAbove * exp ( -IM*AAboveView[localIndex]*h ) ) * IM/ ( sqrt ( 2 ) *2*h );

        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONCAVE:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiKLeft  * exp ( IM*ALeftView[localIndex] *h )
                + psiKRight * exp ( -IM*ARightView[localIndex]*h )
                - psiKBelow * exp ( IM*ABelowView[localIndex]*h )
                + psiKAbove * exp ( -IM*AAboveView[localIndex]*h ) ) * IM/ ( sqrt ( 2 ) *2*h );

        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONCAVE:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( + psiKLeft  * exp ( IM*ALeftView[localIndex] *h )
                - psiKRight * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow * exp ( IM*ABelowView[localIndex]*h )
                - psiKAbove * exp ( -IM*AAboveView[localIndex]*h ) ) * IM/ ( sqrt ( 2 ) *2*h );

        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONCAVE:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiKLeft  * exp ( IM*ALeftView[localIndex] *h )
                + psiKRight * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow * exp ( IM*ABelowView[localIndex]*h )
                - psiKAbove * exp ( -IM*AAboveView[localIndex]*h ) ) * IM/ ( sqrt ( 2 ) *2*h );

        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;
    default:
        TEST_FOR_EXCEPTION ( true,
#include "Ginla_Komplex_DoubleMatrix.h"
                             std::logic_error,
                             "Illegal not type \"" << nt << "\"." );
    }

    // return the result
    return res;
}
// =============================================================================
Teuchos::RCP<const Ginla::State::Virtual>
Ginla::FDM::Operator::BCCentral::
getDFDh0( const Teuchos::RCP<const Ginla::FDM::State> & state ) const
{ 
  // initialize F
  Teuchos::RCP<Ginla::State::Virtual> F =
      Teuchos::rcp( new Ginla::FDM::State( state->getPsi()->getMap(),
                                           state->getGrid() ) );
                                      
  Teuchos::ArrayRCP<double_complex> FView = F->getPsiNonConst()->get1dViewNonConst();

  // rebuild cache for A?
  if ( cacheNeedsUpdating_ )
    this->buildACache_();
  
  Teuchos::ArrayRCP<const double> ALeftView  = ALeft_->get1dView();
  Teuchos::ArrayRCP<const double> ARightView = ARight_->get1dView();
  Teuchos::ArrayRCP<const double> ABelowView = ABelow_->get1dView();
  Teuchos::ArrayRCP<const double> AAboveView = AAbove_->get1dView();
  Teuchos::ArrayRCP<const double> dAdH0LeftView  = dAdH0Left_->get1dView();
  Teuchos::ArrayRCP<const double> dAdH0RightView = dAdH0Right_->get1dView();
  Teuchos::ArrayRCP<const double> dAdH0BelowView = dAdH0Below_->get1dView();
  Teuchos::ArrayRCP<const double> dAdH0AboveView = dAdH0Above_->get1dView();
  
  Teuchos::ArrayRCP<const double_complex> psiView = state->getPsi()->get1dView();
  double chi = state->getChi();
  
  double h = grid_->getUniformH();
  
  // loop over the nodes
  unsigned int localLength = F->getPsi()->getLocalLength();
  
  for ( unsigned int localIndex=0; localIndex<localLength; localIndex++ )
  {
      int globalIndex = rangeMap_->getGlobalElement ( localIndex );
      
      FView[localIndex] = this->getDFDh0Entry_( psiView,
                                                chi,
                                                localIndex,
                                                globalIndex,
                                                ALeftView,
                                                ARightView,
                                                ABelowView,
                                                AAboveView,
                                                dAdH0LeftView,
                                                dAdH0RightView,
                                                dAdH0BelowView,
                                                dAdH0AboveView,
                                                h
                                              );
  }
  
  return F;
}
// =============================================================================
double_complex
Ginla::FDM::Operator::BCCentral::
getDFDh0Entry_( Teuchos::ArrayRCP<const double_complex> & psiView,
                const double                              chi,
                const int localIndex,
                const int globalIndex,
                Teuchos::ArrayRCP<const double> & ALeftView,
                Teuchos::ArrayRCP<const double> & ARightView,
                Teuchos::ArrayRCP<const double> & ABelowView,
                Teuchos::ArrayRCP<const double> & AAboveView,
                Teuchos::ArrayRCP<const double> & dAdH0LeftView,
                Teuchos::ArrayRCP<const double> & dAdH0RightView,
                Teuchos::ArrayRCP<const double> & dAdH0BelowView,
                Teuchos::ArrayRCP<const double> & dAdH0AboveView,
                const double                      h
              ) const
{
    double_complex res;
    double_complex psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;

    Recti::Grid::Abstract::nodeType nt = grid_->getNodeType ( globalIndex );    
    
    psiK = psiView[localIndex];
    
    switch ( nt )
    {
    case Recti::Grid::Abstract::INTERIOR:
        // TODO Gets the local index. ==> Only works on one core.
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( + psiKLeft * IM*dAdH0LeftView[localIndex] *h * exp ( IM*ALeftView[localIndex] *h ) - psiKRight* IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow* IM*dAdH0BelowView[localIndex]*h * exp ( IM*ABelowView[localIndex]*h ) - psiKAbove* IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) / ( h*h );

        res *= exp ( IM*chi );
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONVEX:
        // -------------------------------------------------------------------
        // interior equation, then outward derivative substituted
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( - psiKRight* 2.0 * IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                - psiKAbove* 2.0 * IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) / ( h*h );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONVEX:
        // -----------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( + psiKLeft * 2.0 * IM*dAdH0LeftView[localIndex] *h * exp ( IM*ALeftView[localIndex] *h )
                - psiKAbove* 2.0 * IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) / ( h*h );
        res *= exp ( IM*chi );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONVEX:
        // -----------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( + psiKLeft * 2.0 * IM*dAdH0LeftView[localIndex] *h * exp ( IM*ALeftView[localIndex] *h ) 
                + psiKBelow* 2.0 * IM*dAdH0BelowView[localIndex]*h * exp ( IM*ABelowView[localIndex]*h )
              ) / ( h*h );
        res *= exp ( IM*chi );
        // -----------------------------------------------------------------------

        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONVEX:
        // -----------------------------------------------------------------------
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( - psiKRight* 2.0 * IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow* 2.0 * IM*dAdH0BelowView[localIndex]*h * exp (  IM*ABelowView[localIndex]*h )
              ) / ( h*h );
        res *= exp ( IM*chi );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOM:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( + psiKLeft *       IM*dAdH0LeftView[localIndex] *h * exp (  IM*ALeftView[localIndex] *h ) - psiKRight       * IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                - psiKAbove* 2.0 * IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) / ( h*h );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_RIGHT:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( + psiKLeft * 2.0 * IM*dAdH0LeftView[localIndex] *h * exp ( IM*ALeftView[localIndex] *h )
                + psiKBelow      * IM*dAdH0BelowView[localIndex]*h * exp ( IM*ABelowView[localIndex]*h ) - psiKAbove * IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) / ( h*h );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOP:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( + psiKLeft *       IM*dAdH0LeftView[localIndex] *h * exp ( IM*ALeftView[localIndex] *h ) - psiKRight* IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow* 2.0 * IM*dAdH0BelowView[localIndex]*h * exp ( IM*ABelowView[localIndex]*h )
              ) / ( h*h );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_LEFT:
        // -------------------------------------------------------------------
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( - psiKRight* 2.0 * IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow*       IM*dAdH0BelowView[localIndex]*h * exp (  IM*ABelowView[localIndex]*h ) - psiKAbove* IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) / ( h*h );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONCAVE:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( + psiKLeft  * IM*dAdH0LeftView[localIndex] *h * exp ( IM*ALeftView[localIndex] *h )
                + psiKRight * IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                - psiKBelow * IM*dAdH0BelowView[localIndex]*h * exp ( IM*ABelowView[localIndex]*h )
                - psiKAbove * IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) * IM/ ( sqrt ( 2 ) *2*h );

        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONCAVE:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( - psiKLeft  * IM*dAdH0LeftView[localIndex] *h * exp ( IM*ALeftView[localIndex] *h )
                - psiKRight * IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                - psiKBelow * IM*dAdH0BelowView[localIndex]*h * exp ( IM*ABelowView[localIndex]*h )
                - psiKAbove * IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) * IM/ ( sqrt ( 2 ) *2*h );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONCAVE:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( + psiKLeft  * IM*dAdH0LeftView[localIndex] *h * exp ( IM*ALeftView[localIndex] *h )
                + psiKRight * IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow * IM*dAdH0BelowView[localIndex]*h * exp ( IM*ABelowView[localIndex]*h )
                + psiKAbove * IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) * IM/ ( sqrt ( 2 ) *2*h );

        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONCAVE:
        // -------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        res = ( - psiKLeft  * IM*dAdH0LeftView[localIndex] *h * exp ( IM*ALeftView[localIndex] *h )
                - psiKRight * IM*dAdH0RightView[localIndex]*h * exp ( -IM*ARightView[localIndex]*h )
                + psiKBelow * IM*dAdH0BelowView[localIndex]*h * exp ( IM*ABelowView[localIndex]*h )
                + psiKAbove * IM*dAdH0AboveView[localIndex]*h * exp ( -IM*AAboveView[localIndex]*h )
              ) * IM/ ( sqrt ( 2 ) *2*h );

        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;
    default:
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Illegal not type \"" << nt << "\"." );
    }

    // return the result
    return res;
}
// =============================================================================
Teuchos::RCP<const Ginla::Komplex::DoubleMatrix>
Ginla::FDM::Operator::BCCentral::
getJacobian ( const Teuchos::RCP<const Ginla::FDM::State> & state
            )
{
  if ( firstTime_ )
  {
      const Teuchos::RCP<const ComplexMap> map = state->getPsi()->getMap();
      AB_ = Teuchos::rcp( new Ginla::Komplex::DoubleMatrix( map, map ) );
  }

  // rebuild cache for A?
  if ( cacheNeedsUpdating_ )
    this->buildACache_();
  
  Teuchos::ArrayRCP<const double> ALeftView  = ALeft_->get1dView();
  Teuchos::ArrayRCP<const double> ARightView = ARight_->get1dView();
  Teuchos::ArrayRCP<const double> ABelowView = ABelow_->get1dView();
  Teuchos::ArrayRCP<const double> AAboveView = AAbove_->get1dView();

  Teuchos::ArrayRCP<const double_complex> psiView = state->getPsi()->get1dView();
  double chi = state->getChi();

  double h = grid_->getUniformH();

  // loop over the nodes
  unsigned int localLength = state->getPsi()->getLocalLength();
  
  Teuchos::Array<Thyra::Ordinal> columnIndicesPsi;
  Teuchos::Array<double_complex> valuesPsi;
  Teuchos::Array<Thyra::Ordinal> columnIndicesPsiConj;
  Teuchos::Array<double_complex> valuesPsiConj;

  for ( unsigned int localRow=0; localRow<localLength; localRow++ )
  {
      int globalRow= rangeMap_->getGlobalElement ( localRow );

      this->getJacobianRow_ ( psiView,
                              chi,
                              localRow,
                              globalRow,
                              ALeftView, ARightView, ABelowView, AAboveView,
                              h,
                              columnIndicesPsi, valuesPsi,
                              columnIndicesPsiConj, valuesPsiConj
                            );

      // add the rows to the Jacobian
      AB_->putALocalValues( localRow,
                            columnIndicesPsi(),
                            valuesPsi() );

      AB_->putBLocalValues( localRow,
                            columnIndicesPsiConj(),
                            valuesPsiConj() );
  }

  if ( firstTime_ )
  {
    AB_->finalize();
    firstTime_ = false;
  }

  return AB_;
}
// =============================================================================
void
Ginla::FDM::Operator::BCCentral::
getJacobianRow_ ( Teuchos::ArrayRCP<const double_complex> & psiView,
                  const double                              chi,
                  const int localIndex,
                  const int globalIndex,
                  Teuchos::ArrayRCP<const double> & ALeftView,
                  Teuchos::ArrayRCP<const double> & ARightView,
                  Teuchos::ArrayRCP<const double> & ABelowView,
                  Teuchos::ArrayRCP<const double> & AAboveView,
                  const double                      h,
                  Teuchos::Array<Thyra::Ordinal>  & columnIndicesPsi,
                  Teuchos::Array<double_complex>  & valuesPsi,
                  Teuchos::Array<Thyra::Ordinal>  & columnIndicesPsiConj,
                  Teuchos::Array<double_complex>  & valuesPsiConj
                ) const
{
    int kLeft, kRight, kBelow, kAbove;
    int numEntriesPsi, numEntriesPsiConj;
    
    Recti::Grid::Abstract::nodeType nt = grid_->getNodeType ( globalIndex );
    
    switch ( nt )
    {
    case Recti::Grid::Abstract::INTERIOR:
        // ---------------------------------------------------------------------
        kRight = grid_->getKRight ( globalIndex );
        kLeft  = grid_->getKLeft ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );
        kBelow = grid_->getKBelow ( globalIndex );

        numEntriesPsi = 5;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = globalIndex;
        columnIndicesPsi[1] = kLeft;
        columnIndicesPsi[2] = kRight;
        columnIndicesPsi[3] = kBelow;
        columnIndicesPsi[4] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = - 4.0            / ( h*h )
                       + ( 1 - 2.0*norm ( psiView[localIndex] ) );
        valuesPsi[1] = exp ( IM*ALeftView[localIndex] *h ) / ( h*h );
        valuesPsi[2] = exp ( -IM*ARightView[localIndex]*h ) / ( h*h );
        valuesPsi[3] = exp ( IM*ABelowView[localIndex]*h ) / ( h*h );
        valuesPsi[4] = exp ( -IM*AAboveView[localIndex]*h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = globalIndex;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[localIndex]*psiView[localIndex];

        break;
    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONVEX:
        // -------------------------------------------------------------------
        kRight = grid_->getKRight ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );

        numEntriesPsi = 3;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = globalIndex;
        columnIndicesPsi[1] = kRight;
        columnIndicesPsi[2] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = - 4.0 / ( h*h )
                       + ( 1 - 2.0*norm ( psiView[localIndex] ) );
        valuesPsi[1] = 2.0 * exp ( -IM*ARightView[localIndex]*h ) / ( h*h );
        valuesPsi[2] = 2.0 * exp ( -IM*AAboveView[localIndex]*h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = globalIndex;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[localIndex]*psiView[localIndex];
        valuesPsiConj[0] *= exp ( IM*chi*2.0 );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONVEX:
        // -----------------------------------------------------------------------
        kLeft  = grid_->getKLeft ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );

        numEntriesPsi = 3;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = globalIndex;
        columnIndicesPsi[1] = kLeft;
        columnIndicesPsi[2] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -4.0 / ( h*h )
                       + ( 1 - 2.0*norm ( psiView[localIndex] ) );
        valuesPsi[1] = 2.0 * exp ( IM*ALeftView[localIndex] *h ) / ( h*h );
        valuesPsi[2] = 2.0 * exp ( -IM*AAboveView[localIndex]*h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = globalIndex;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[localIndex]*psiView[localIndex];
        valuesPsiConj[0] *= exp ( IM*chi*2.0 );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONVEX:
        // -----------------------------------------------------------------------
        kLeft  = grid_->getKLeft ( globalIndex );
        kBelow = grid_->getKBelow ( globalIndex );

        numEntriesPsi = 3;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = globalIndex;
        columnIndicesPsi[1] = kLeft;
        columnIndicesPsi[2] = kBelow;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -4.0 / ( h*h )
                       + ( 1 - 2.0*norm ( psiView[localIndex] ) );
        valuesPsi[1] = 2.0 * exp ( IM*ALeftView[localIndex] *h ) / ( h*h );
        valuesPsi[2] = 2.0 * exp ( IM*ABelowView[localIndex]*h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = globalIndex;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[localIndex]*psiView[localIndex];
        valuesPsiConj[0] *= exp ( IM*chi*2.0 );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONVEX:
        // -----------------------------------------------------------------------
        kRight = grid_->getKRight ( globalIndex );
        kBelow = grid_->getKBelow ( globalIndex );

        numEntriesPsi = 3;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = globalIndex;
        columnIndicesPsi[1] = kRight;
        columnIndicesPsi[2] = kBelow;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -4.0 / ( h*h )
                       + ( 1 - 2.0*norm ( psiView[localIndex] ) );
        valuesPsi[1] = 2.0 * exp ( -IM*ARightView[localIndex]*h ) / ( h*h );
        valuesPsi[2] = 2.0 * exp ( IM*ABelowView[localIndex]*h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = globalIndex;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[localIndex]*psiView[localIndex];
        valuesPsiConj[0] *= exp ( IM*chi*2.0 );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOM:
        // -------------------------------------------------------------------
        kLeft  = grid_->getKLeft ( globalIndex );
        kRight = grid_->getKRight ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );

        numEntriesPsi = 4;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = globalIndex;
        columnIndicesPsi[1] = kLeft;
        columnIndicesPsi[2] = kRight;
        columnIndicesPsi[3] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = - 4.0            / ( h*h )
                       + ( 1.0 - 2.0*norm ( psiView[localIndex] ) );
        valuesPsi[1] =       exp ( IM*ALeftView[localIndex] *h ) / ( h*h );
        valuesPsi[2] =       exp ( -IM*ARightView[localIndex]*h ) / ( h*h );
        valuesPsi[3] = 2.0 * exp ( -IM*AAboveView[localIndex]*h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = globalIndex;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[localIndex]*psiView[localIndex];
        valuesPsiConj[0] *= exp ( IM*chi*2.0 );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_RIGHT:
        // -------------------------------------------------------------------
        kBelow = grid_->getKBelow ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );
        kLeft  = grid_->getKLeft ( globalIndex );

        numEntriesPsi = 4;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = globalIndex;
        columnIndicesPsi[1] = kBelow;
        columnIndicesPsi[2] = kAbove;
        columnIndicesPsi[3] = kLeft;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = - 4.0            / ( h*h )
                       + ( 1.0 - 2.0*norm ( psiView[localIndex] ) );
        valuesPsi[1] =       exp ( IM*ABelowView[localIndex]*h ) / ( h*h );
        valuesPsi[2] =       exp ( -IM*AAboveView[localIndex]*h ) / ( h*h );
        valuesPsi[3] = 2.0 * exp ( IM*ALeftView[localIndex] *h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = globalIndex;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[localIndex]*psiView[localIndex];
        valuesPsiConj[0] *= exp ( IM*chi*2.0 );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOP:
        // -------------------------------------------------------------------
        kBelow = grid_->getKBelow ( globalIndex );
        kRight = grid_->getKRight ( globalIndex );
        kLeft  = grid_->getKLeft ( globalIndex );

        numEntriesPsi = 4;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = globalIndex;
        columnIndicesPsi[1] = kBelow;
        columnIndicesPsi[2] = kLeft;
        columnIndicesPsi[3] = kRight;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = - 4.0             / ( h*h )
                       + ( 1.0 - 2.0*norm ( psiView[localIndex] ) );
        valuesPsi[1] = 2.0 * exp ( IM*ABelowView[localIndex]*h ) / ( h*h );
        valuesPsi[2] =       exp ( IM*ALeftView[localIndex] *h ) / ( h*h );
        valuesPsi[3] =       exp ( -IM*ARightView[localIndex]*h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = globalIndex;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[localIndex]*psiView[localIndex];
        valuesPsiConj[0] *= exp ( IM*chi*2.0 );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_LEFT:
        // -------------------------------------------------------------------
        kBelow = grid_->getKBelow ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );
        kRight = grid_->getKRight ( globalIndex );

        numEntriesPsi = 4;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = globalIndex;
        columnIndicesPsi[1] = kBelow;
        columnIndicesPsi[2] = kAbove;
        columnIndicesPsi[3] = kRight;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = - 4.0            / ( h*h )
                       + ( 1.0 - 2.0*norm ( psiView[localIndex] ) );
        valuesPsi[1] =       exp ( IM*ABelowView[localIndex]*h ) / ( h*h );
        valuesPsi[2] =       exp ( -IM*AAboveView[localIndex]*h ) / ( h*h );
        valuesPsi[3] = 2.0 * exp ( -IM*ARightView[localIndex]*h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = globalIndex;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[localIndex]*psiView[localIndex];
        valuesPsiConj[0] *= exp ( IM*chi*2.0 );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONCAVE:
        // -----------------------------------------------------------------------
        kRight = grid_->getKRight ( globalIndex );
        kLeft  = grid_->getKLeft ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );
        kBelow = grid_->getKBelow ( globalIndex );

        numEntriesPsi = 4;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = kLeft;
        columnIndicesPsi[1] = kRight;
        columnIndicesPsi[2] = kBelow;
        columnIndicesPsi[3] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] =  exp ( IM*ALeftView[localIndex]  *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[1] = -exp ( -IM*ARightView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[2] = -exp ( IM*ABelowView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[3] =  exp ( -IM*AAboveView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONCAVE:
        // -----------------------------------------------------------------------
        kRight = grid_->getKRight ( globalIndex );
        kLeft  = grid_->getKLeft ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );
        kBelow = grid_->getKBelow ( globalIndex );

        numEntriesPsi = 4;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = kLeft;
        columnIndicesPsi[1] = kRight;
        columnIndicesPsi[2] = kBelow;
        columnIndicesPsi[3] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -exp ( IM*ALeftView[localIndex]  *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[1] =  exp ( -IM*ARightView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[2] = -exp ( IM*ABelowView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[3] =  exp ( -IM*AAboveView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONCAVE:
        // -----------------------------------------------------------------------
        kRight = grid_->getKRight ( globalIndex );
        kLeft  = grid_->getKLeft ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );
        kBelow = grid_->getKBelow ( globalIndex );

        numEntriesPsi = 4;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = kLeft;
        columnIndicesPsi[1] = kRight;
        columnIndicesPsi[2] = kBelow;
        columnIndicesPsi[3] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] =  exp ( IM*ALeftView[localIndex]  *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[1] = -exp ( -IM*ARightView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[2] =  exp ( IM*ABelowView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[3] = -exp ( -IM*AAboveView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONCAVE:
        // -----------------------------------------------------------------------
        kRight = grid_->getKRight ( globalIndex );
        kLeft  = grid_->getKLeft ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );
        kBelow = grid_->getKBelow ( globalIndex );

        numEntriesPsi = 4;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = kLeft;
        columnIndicesPsi[1] = kRight;
        columnIndicesPsi[2] = kBelow;
        columnIndicesPsi[3] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -exp ( IM*ALeftView[localIndex]  *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[1] =  exp ( -IM*ARightView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[2] =  exp ( IM*ABelowView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[3] = -exp ( -IM*AAboveView[localIndex] *h ) * IM/ ( sqrt ( 2 ) *2*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
        // -----------------------------------------------------------------------
        break;

    default:
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Illegal note type \"" << nt << "\"." );
    }

}
// =============================================================================
