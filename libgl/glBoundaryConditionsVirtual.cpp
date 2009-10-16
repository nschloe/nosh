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
void
GlBoundaryConditionsVirtual::getEquationType ( const int eqIndex,
                                               const StaggeredGrid::StaggeredGrid &sGrid,
                                               GlBoundaryConditionsVirtual::equationType  &eqType,
                                               Teuchos::Array<int>                &i       )
{
  int Nx = sGrid.getNx();

  if ( eqIndex==0 )
    {
      eqType = GlBoundaryConditionsVirtual::BOTTOMLEFT;
      i[0] = 0;
      i[1] = 0;
    }
  else if ( eqIndex==Nx )
    {
      eqType = GlBoundaryConditionsVirtual::BOTTOMRIGHT;
      i[0] = Nx;
      i[1] = 0;
    }
  else if ( eqIndex==2*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::TOPRIGHT;
      i[0] = Nx;
      i[1] = Nx;
    }
  else if ( eqIndex==3*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::TOPLEFT;
      i[0] = 0;
      i[1] = Nx;
    }
  else if ( eqIndex<Nx )
    {
      eqType = GlBoundaryConditionsVirtual::BOTTOM;
      i[0] = eqIndex;
      i[1] = 0;
    }
  else if ( eqIndex<2*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::RIGHT;
      i[0] = Nx;
      i[1] = eqIndex-Nx;
    }
  else if ( eqIndex<3*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::TOP;
      i[0] = 3*Nx-eqIndex;
      i[1] = Nx;
    }
  else if ( eqIndex<4*Nx )
    {
      eqType = GlBoundaryConditionsVirtual::LEFT;
      i[0] = 0;
      i[1] = 4*Nx-eqIndex;
    }
  else
    {
      throw glException ( "GlBoundaryConditionsVirtual::getEquationType",
                          "Illegal running index eqIndex="
                          + EpetraExt::toString ( eqIndex ) );
    }
}
// =============================================================================