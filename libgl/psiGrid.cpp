#include "psiGrid.h"

// =============================================================================
// Class constructor
PsiGrid::PsiGrid( int nx ):
Nx(nx),
d(2)
{
}
// =============================================================================


// =============================================================================
// Destructor
PsiGrid::~PsiGrid()
{
}
// =============================================================================


// =============================================================================
// maps a running index k to a 2D index i
int* PsiGrid::k2i( int k )
{
  static int i[2];

  if (k<Nx) { // lower shore
    i[0] = k-1;
    i[1] = 0;
  }
  else if (k<2*Nx) { // right shore
    i[0] = Nx;
    i[1] = k-Nx;
  }
  else if (k<3*Nx) { // upper shore
    i[0] = 3*Nx-k;
    i[1] = Nx;
  }
  else if (k<4*Nx) { // left shore
    i[0] = 0;
    i[1] = 4*Nx-k;
  }
  else { // on the interior
    int numBoundaryNodes = 4*Nx;
    i[0] = (k-numBoundaryNodes)%Nx + 1;
    i[1] = (k-numBoundaryNodes) / Nx;
  }

  return i;
}
// =============================================================================

