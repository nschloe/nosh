#include "glBoundaryConditionsVirtual.h"
#include "glException.h"

#include <EpetraExt_Utils.h>

// =============================================================================
GlBoundaryConditionsVirtual::GlBoundaryConditionsVirtual()
{
}
// =============================================================================
GlBoundaryConditionsVirtual::~GlBoundaryConditionsVirtual()
{
}
// =============================================================================
//! With the \cboundaryEquationIndex-th boundary equation, this functions
//! connects the equations ``centered'' around the \cboundaryEquationIndex-th
//! boundary node.
void
GlBoundaryConditionsVirtual::getEquationType ( const int boundaryEquationIndex,
                                               const Grid::Grid  &grid,
                                               GlBoundaryConditionsVirtual::equationType  &eqType,
                                               Teuchos::Array<int>                &i       )
{

//  Grid::notType nt = grid.boundaryNodeType( boundaryEquationIndex );
// the switch to translate
// also return index K around with it is centered!

  int Nx = grid.getNx();

  if ( boundaryEquationIndex==0 )
    {
      eqType = GlBoundaryConditionsVirtual::BOTTOMLEFT;
      i[0] = 0;
      i[1] = 0;
    }
  else if ( boundaryEquationIndex==Nx )
    {
      eqType = GlBoundaryConditionsVirtual::BOTTOMRIGHT;
      i[0] = Nx;
      i[1] = 0;
    }
  else if ( boundaryEquationIndex==2*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::TOPRIGHT;
      i[0] = Nx;
      i[1] = Nx;
    }
  else if ( boundaryEquationIndex==3*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::TOPLEFT;
      i[0] = 0;
      i[1] = Nx;
    }
  else if ( boundaryEquationIndex<Nx )
    {
      eqType = GlBoundaryConditionsVirtual::BOTTOM;
      i[0] = boundaryEquationIndex;
      i[1] = 0;
    }
  else if ( boundaryEquationIndex<2*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::RIGHT;
      i[0] = Nx;
      i[1] = boundaryEquationIndex-Nx;
    }
  else if ( boundaryEquationIndex<3*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::TOP;
      i[0] = 3*Nx-boundaryEquationIndex;
      i[1] = Nx;
    }
  else if ( boundaryEquationIndex<4*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::LEFT;
      i[0] = 0;
      i[1] = 4*Nx-boundaryEquationIndex;
    }
  else
    {
      throw glException ( "GlBoundaryConditionsVirtual::getEquationType",
                          "Illegal running index boundaryEquationIndex="
                          + EpetraExt::toString ( boundaryEquationIndex ) );
    }
}
// =============================================================================
