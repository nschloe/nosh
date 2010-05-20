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

#include "Ginla_Operator_BCInner.h"

#include "Recti_Grid_Uniform.h"
#include "Ginla_MagneticVectorPotential_Centered.h"

#include <Teuchos_Array.hpp>


// =============================================================================
Ginla::Operator::BCInner::
BCInner ( const Teuchos::RCP<Recti::Grid::Uniform>                     & grid,
          const Teuchos::RCP<Ginla::MagneticVectorPotential::Centered> & A,
          const Teuchos::RCP<const ComplexMap>                         & domainMap,
          const Teuchos::RCP<const ComplexMap>                         & rangeMap
        ) :
        Ginla::Operator::Virtual ( grid, A, domainMap, rangeMap )
{
}
// =============================================================================
Ginla::Operator::BCInner::~BCInner()
{
}
// =============================================================================
Teuchos::RCP<Ginla::State>
Ginla::Operator::BCInner::
getF( const Teuchos::RCP<const Ginla::State> & state ) const
{
  // initialize F
  Teuchos::RCP<Ginla::State> F =
      Teuchos::rcp( new Ginla::State( state->getPsi()->getMap(),
                                      state->getGrid() ) );
                                      
  Teuchos::ArrayRCP<double_complex> FView = F->getPsiNonConst()->get1dViewNonConst();
  
  // loop over the nodes
  unsigned int localLength = F->getPsi()->getLocalLength();
  
  for ( unsigned int k=0; k<localLength; k++ )
  {
      int globalIndex = rangeMap_->getGlobalElement ( k );
      FView[k] = this->getFEntry_ ( state, k, globalIndex );
  }
  
  return F;
}
// =============================================================================
double_complex
Ginla::Operator::BCInner::
getFEntry_ ( const Teuchos::RCP<const Ginla::State> & state,
             const int localIndex,
             const int globalIndex
           ) const
{
    double_complex res;
    double_complex psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;
    double ARightView[k], ALeftView[k], AAboveView[k], ABelowView[k];

    double h = grid_->getUniformH();

    Recti::Grid::Abstract::nodeType nt = grid_->getNodeType ( k );

    // Get a view of the whole vector.
    // Remember: This only works with one core.
    Teuchos::ArrayRCP<const double_complex> psiView = state->getPsi()->get1dView();
    double chi = state->getChi();

    psiK = psiView[k];
    
    Teuchos::ArrayRCP<double> ALeftView  = ALeft_->get1dViewNonConst();
    Teuchos::ArrayRCP<double> ARightView = ARight_->get1dViewNonConst();
    Teuchos::ArrayRCP<double> ABelowView = ABelow_->get1dViewNonConst();
    Teuchos::ArrayRCP<double> AAboveView = AAbove_->get1dViewNonConst();

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
                + psiKLeft*  exp ( IM*ALeftView[k] *h ) + psiKRight* exp ( -IM*ARightView[k]*h )
                + psiKBelow* exp ( IM*ABelowView[k]*h ) + psiKAbove* exp ( -IM*AAboveView[k]*h ) ) / ( h*h )
              + psiK * ( 1-norm ( psiK ) );
        res *= exp ( IM*chi );

        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONVEX:
    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONCAVE:
        // -------------------------------------------------------------------
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiK      * 2.0
                + psiKRight * exp ( -IM*ARightView[k]*h )
                + psiKAbove * exp ( -IM*AAboveView[k]*h ) ) * IM/ ( sqrt ( 2 ) *h );
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONVEX:
    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONCAVE:
        // -----------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiK      * 2.0
                + psiKLeft  * exp ( IM*ALeftView[k] *h )
                + psiKAbove * exp ( -IM*AAboveView[k]*h ) ) * IM/ ( sqrt ( 2 ) *h );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONVEX:
    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONCAVE:
        // -----------------------------------------------------------------------
        psiKLeft  = psiView[ grid_->getKLeft ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];

        res = ( - psiK      * 2.0
                + psiKLeft  * exp ( IM*ALeftView[k] *h )
                + psiKBelow * exp ( IM*ABelowView[k]*h ) ) * IM/ ( sqrt ( 2 ) *h );
        res *= exp ( IM*chi );
        // -----------------------------------------------------------------------

        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONVEX:
    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONCAVE:
        // -----------------------------------------------------------------------
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];

        res = ( - psiK      * 2.0
                + psiKRight * exp ( -IM*ARightView[k]*h )
                + psiKBelow * exp ( IM*ABelowView[k]*h ) ) * IM/ ( sqrt ( 2 ) *h );
        res *= exp ( IM*chi );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOM:
        // -------------------------------------------------------------------
        // normal derivative
        psiKAbove = psiView[ grid_->getKAbove ( globalIndex ) ];

        res = ( - psiK
                + psiKAbove * exp ( -IM*AAboveView[k]*h ) ) * IM/h;
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_RIGHT:
        // -------------------------------------------------------------------
        // normal derivative
        psiKLeft = psiView[ grid_->getKLeft ( globalIndex ) ];

        res = ( - psiK
                + psiKLeft * exp ( IM*ALeftView[k]*h ) ) * IM/h;
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOP:
        // -------------------------------------------------------------------
        // normal derivative
        psiKBelow = psiView[ grid_->getKBelow ( globalIndex ) ];

        res = ( - psiK
                + psiKBelow * exp ( IM*ABelowView[k]*h ) ) * IM/h;
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_LEFT:
        // -------------------------------------------------------------------
        // normal derivative
        psiKRight = psiView[ grid_->getKRight ( globalIndex ) ];

        res = ( - psiK
                + psiKRight * exp ( -IM*ARightView[k]*h ) ) * IM/h;
        res *= exp ( IM*chi );
        // -------------------------------------------------------------------
        break;

    default:
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Illegal node type \"" << nt << "\"." );
    }

    // return the result
    return res;
}
// =============================================================================
void
Ginla::Operator::BCInner::
getJacobianRow ( const Teuchos::RCP<const Ginla::State> & state,
                 const int                                k,
                 const int                                globalIndex,
                 Teuchos::Array<int>                    & columnIndicesPsi,
                 Teuchos::Array<double_complex>         & valuesPsi,
                 Teuchos::Array<int>                    & columnIndicesPsiConj,
                 Teuchos::Array<double_complex>         & valuesPsiConj
               ) const
{
    int kLeft, kRight, kBelow, kAbove;
    int numEntriesPsi, numEntriesPsiConj;

    double h = grid_->getUniformH();
    
    Teuchos::ArrayRCP<double> ALeftView  = ALeft_->get1dViewNonConst();
    Teuchos::ArrayRCP<double> ARightView = ARight_->get1dViewNonConst();
    Teuchos::ArrayRCP<double> ABelowView = ABelow_->get1dViewNonConst();
    Teuchos::ArrayRCP<double> AAboveView = AAbove_->get1dViewNonConst();

    Teuchos::ArrayRCP<const double_complex> psiView = state->getPsi()->get1dView();

    Recti::Grid::Abstract::nodeType nt = grid_->getNodeType(k);
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
        columnIndicesPsi[0] = k;
        columnIndicesPsi[1] = kLeft;
        columnIndicesPsi[2] = kRight;
        columnIndicesPsi[3] = kBelow;
        columnIndicesPsi[4] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = - 4.0            / ( h*h )
                       + ( 1 - 2.0*norm ( psiView[k] ) );
        valuesPsi[1] = exp ( IM*ALeftView[k] *h ) / ( h*h );
        valuesPsi[2] = exp ( -IM*ARightView[k]*h ) / ( h*h );
        valuesPsi[3] = exp ( IM*ABelowView[k]*h ) / ( h*h );
        valuesPsi[4] = exp ( -IM*AAboveView[k]*h ) / ( h*h );

        numEntriesPsiConj = 1;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        columnIndicesPsiConj[0] = k;

        valuesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj[0] = -psiView[k]*psiView[k];

        break;
    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONVEX:
    case Recti::Grid::Abstract::BOUNDARY_BOTTOMLEFTCONCAVE:
        // -------------------------------------------------------------------
        kRight = grid_->getKRight ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );

        numEntriesPsi = 3;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = k;
        columnIndicesPsi[1] = kRight;
        columnIndicesPsi[2] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -2.0                * IM/ ( sqrt ( 2 ) *h );
        valuesPsi[1] = exp ( -IM*ARightView[k]*h ) * IM/ ( sqrt ( 2 ) *h );
        valuesPsi[2] = exp ( -IM*AAboveView[k]*h ) * IM/ ( sqrt ( 2 ) *h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONVEX:
    case Recti::Grid::Abstract::BOUNDARY_BOTTOMRIGHTCONCAVE:
        // -----------------------------------------------------------------------
        kLeft  = grid_->getKLeft ( globalIndex );
        kAbove = grid_->getKAbove ( globalIndex );

        numEntriesPsi = 3;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = k;
        columnIndicesPsi[1] = kLeft;
        columnIndicesPsi[2] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -2.0                * IM/ ( sqrt ( 2 ) *h );
        valuesPsi[1] = exp ( IM*ALeftView[k] *h )  * IM/ ( sqrt ( 2 ) *h );
        valuesPsi[2] = exp ( -IM*AAboveView[k]*h ) * IM/ ( sqrt ( 2 ) *h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );

        valuesPsiConj.resize ( numEntriesPsiConj );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONVEX:
    case Recti::Grid::Abstract::BOUNDARY_TOPRIGHTCONCAVE:
        // -----------------------------------------------------------------------
        kLeft  = grid_->getKLeft ( k );
        kBelow = grid_->getKBelow ( k );

        numEntriesPsi = 3;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = k;
        columnIndicesPsi[1] = kLeft;
        columnIndicesPsi[2] = kBelow;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -2.0             * IM/ ( sqrt ( 2 ) *h );
        valuesPsi[1] = exp ( IM*ALeftView[k] *h ) * IM/ ( sqrt ( 2 ) *h );
        valuesPsi[2] = exp ( IM*ABelowView[k]*h ) * IM/ ( sqrt ( 2 ) *h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );

        valuesPsiConj.resize ( numEntriesPsiConj );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONVEX:
    case Recti::Grid::Abstract::BOUNDARY_TOPLEFTCONCAVE:
        // -----------------------------------------------------------------------
        kRight = grid_->getKRight ( globalIndex );
        kBelow = grid_->getKBelow ( globalIndex );

        numEntriesPsi = 3;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = k;
        columnIndicesPsi[1] = kRight;
        columnIndicesPsi[2] = kBelow;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -2.0                * IM/ ( sqrt ( 2 ) *h );
        valuesPsi[1] = exp ( -IM*ARightView[k]*h ) * IM/ ( sqrt ( 2 ) *h );
        valuesPsi[2] = exp ( IM*ABelowView[k]*h ) * IM/ ( sqrt ( 2 ) *h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
        // -----------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_BOTTOM:
        // -------------------------------------------------------------------
        // normal derivative
        kAbove = grid_->getKAbove ( k );

        numEntriesPsi = 2;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = k;
        columnIndicesPsi[1] = kAbove;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -1.0                * IM/h;
        valuesPsi[1] = exp ( -IM*AAboveView[k]*h ) * IM/h;

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );

        valuesPsiConj.resize ( numEntriesPsiConj );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_RIGHT:
        // -------------------------------------------------------------------
        // normal derivative
        kLeft = grid_->getKLeft ( k );

        numEntriesPsi = 2;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = k;
        columnIndicesPsi[1] = kLeft;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -1.0            * IM/h;
        valuesPsi[1] = exp ( IM*ALeftView[k]*h ) * IM/h;

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_TOP:
        // -------------------------------------------------------------------
        // normal derivative
        kBelow = grid_->getKBelow ( k );

        numEntriesPsi = 2;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = k;
        columnIndicesPsi[1] = kBelow;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -1.0             * IM/h;
        valuesPsi[1] = exp ( IM*ABelowView[k]*h ) * IM/h;

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
        // -------------------------------------------------------------------
        break;

    case Recti::Grid::Abstract::BOUNDARY_LEFT:
        // -------------------------------------------------------------------
        // normal derivative
        kRight = grid_->getKRight ( k );

        numEntriesPsi = 2;
        columnIndicesPsi.resize ( numEntriesPsi );
        columnIndicesPsi[0] = k;
        columnIndicesPsi[1] = kRight;

        valuesPsi.resize ( numEntriesPsi );
        valuesPsi[0] = -1.0                * IM/h;
        valuesPsi[1] = exp ( -IM*ARightView[k]*h ) * IM/h;

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );

        valuesPsiConj.resize ( numEntriesPsiConj );
        // -------------------------------------------------------------------
        break;

    default:
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Illegal node type \"" << nt << "\"." );
    }

}
// =============================================================================
