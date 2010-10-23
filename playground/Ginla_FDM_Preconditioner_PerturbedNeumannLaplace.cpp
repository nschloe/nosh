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

#include "Ginla_Preconditioner_PerturbedNeumannLaplace.h"

#include "Recti_Grid_Uniform.h"

// =============================================================================
Ginla::Preconditioner::PerturbedNeumannLaplace::
PerturbedNeumannLaplace ( const Teuchos::RCP<Recti::Grid::Uniform> & grid,
                          const Teuchos::RCP<const ComplexMap>     & domainMap,
                          const Teuchos::RCP<const ComplexMap>     & rangeMap
                        ) :
        Ginla::Preconditioner::Virtual ( grid, domainMap, rangeMap ),
        epsilon( 1.0e-4 )
{
}
// =============================================================================
Ginla::Preconditioner::PerturbedNeumannLaplace::
~PerturbedNeumannLaplace()
{
}
// =============================================================================
Teuchos::RCP<const Komplex2::DoubleMatrix>
Ginla::Preconditioner::PerturbedNeumannLaplace::
getMatrix ( const Teuchos::RCP<const Ginla::State> & state
          )
{ 
  if ( firstTime_ )
  {
      const Teuchos::RCP<const ComplexMap> map = state->getPsi()->getMap();
      AB_ = Teuchos::rcp( new Komplex2::DoubleMatrix( map, map ) );
  }
  
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
      int globalRow = rangeMap_->getGlobalElement ( localRow );

      this->getMatrixRow_ ( localRow,
                            globalRow,
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
Ginla::Preconditioner::PerturbedNeumannLaplace::
getMatrixRow_ ( const int localIndex,
                const int globalIndex,
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
        valuesPsi[0] = - 4.0 / ( h*h );
        valuesPsi[1] = 1.0 / ( h*h );
        valuesPsi[2] = 1.0 / ( h*h );
        valuesPsi[3] = 1.0 / ( h*h );
        valuesPsi[4] = 1.0 / ( h*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );

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
        valuesPsi[0] = - 4.0 / ( h*h );
        valuesPsi[1] = 2.0 / ( h*h );
        valuesPsi[2] = 2.0 / ( h*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
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
        valuesPsi[0] = -4.0 / ( h*h );
        valuesPsi[1] = 2.0 / ( h*h );
        valuesPsi[2] = 2.0 / ( h*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
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
        valuesPsi[0] = -4.0 / ( h*h );
        valuesPsi[1] = 2.0 / ( h*h );
        valuesPsi[2] = 2.0 / ( h*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
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
        valuesPsi[0] = -4.0 / ( h*h );
        valuesPsi[1] = 2.0 / ( h*h );
        valuesPsi[2] = 2.0 / ( h*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
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
        valuesPsi[0] = - 4.0 / ( h*h );
        valuesPsi[1] = 1.0 / ( h*h );
        valuesPsi[2] = 1.0 / ( h*h );
        valuesPsi[3] = 2.0 / ( h*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
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
        valuesPsi[0] = - 4.0            / ( h*h );
        valuesPsi[1] = 1.0 / ( h*h );
        valuesPsi[2] = 1.0 / ( h*h );
        valuesPsi[3] = 2.0 / ( h*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
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
        valuesPsi[0] = - 4.0 / ( h*h );
        valuesPsi[1] = 2.0 / ( h*h );
        valuesPsi[2] = 1.0 / ( h*h );
        valuesPsi[3] = 1.0 / ( h*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
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
        valuesPsi[0] = - 4.0 / ( h*h );
        valuesPsi[1] = 1.0 / ( h*h );
        valuesPsi[2] = 1.0 / ( h*h );
        valuesPsi[3] = 2.0 / ( h*h );

        numEntriesPsiConj = 0;
        columnIndicesPsiConj.resize ( numEntriesPsiConj );
        valuesPsiConj.resize ( numEntriesPsiConj );
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
        valuesPsi[0] =  IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[1] = -IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[2] = -IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[3] =  IM/ ( sqrt ( 2 ) *2*h );

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
        valuesPsi[0] = -IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[1] =  IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[2] = -IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[3] =  IM/ ( sqrt ( 2 ) *2*h );

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
        valuesPsi[0] =  IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[1] = -IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[2] =  IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[3] = -IM/ ( sqrt ( 2 ) *2*h );

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
        valuesPsi[0] = -IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[1] =  IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[2] =  IM/ ( sqrt ( 2 ) *2*h );
        valuesPsi[3] = -IM/ ( sqrt ( 2 ) *2*h );

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
