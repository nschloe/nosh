#include "glBoundaryConditionsVirtual.h"

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
GlBoundaryConditionsVirtual::equationType
GlBoundaryConditionsVirtual::getEquationType ( const int boundaryEquationIndex,
                                               const GridSquare::GridSquare &grid  )
{
  GridVirtual::nodeType nt;
  nt = grid.getBoundaryNodeType( boundaryEquationIndex );

  // translate node type into equation type
  switch (nt)
  {
  case GridVirtual::BOTTOMLEFTCONVEX:
    return GlBoundaryConditionsVirtual::BOTTOMLEFTCONVEX;
  case GridVirtual::BOTTOMLEFTCONCAVE:
    return GlBoundaryConditionsVirtual::BOTTOMLEFTCONCAVE;
  case GridVirtual::TOPLEFTCONVEX:
    return GlBoundaryConditionsVirtual::TOPLEFTCONVEX;
  case GridVirtual::TOPLEFTCONCAVE:
    return GlBoundaryConditionsVirtual::TOPLEFTCONCAVE;
  case GridVirtual::BOTTOMRIGHTCONVEX:
    return GlBoundaryConditionsVirtual::BOTTOMRIGHTCONVEX;
  case GridVirtual::BOTTOMRIGHTCONCAVE:
    return GlBoundaryConditionsVirtual::BOTTOMRIGHTCONCAVE;
  case GridVirtual::TOPRIGHTCONVEX:
    return GlBoundaryConditionsVirtual::TOPRIGHTCONVEX;
  case GridVirtual::TOPRIGHTCONCAVE:
    return GlBoundaryConditionsVirtual::TOPRIGHTCONCAVE;
  case GridVirtual::TOP:
    return GlBoundaryConditionsVirtual::TOP;
  case GridVirtual::BOTTOM:
    return GlBoundaryConditionsVirtual::BOTTOM;
  case GridVirtual::LEFT:
    return GlBoundaryConditionsVirtual::LEFT;
  case GridVirtual::RIGHT:
    return GlBoundaryConditionsVirtual::RIGHT;
  default:
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Illegal nodeType=" << nt << ".");
  }

}
// =============================================================================
