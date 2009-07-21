#include "psiGrid.h"
#include <iostream>

#include <stdlib.h>


// =============================================================================
// Class constructor
PsiGrid::PsiGrid( int nx,
                  double edgelength ):
  Nx(nx),
  Edgelength(edgelength) // not strictly necessary
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
int PsiGrid::getNx()
{
    return Nx;
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
    i[0] = (k-numBoundaryNodes)%(Nx-1) + 1;
    i[1] = (k-numBoundaryNodes)/(Nx-1);
  }

  return i;
}
// =============================================================================


// =============================================================================
// maps a 2D index i to a running index k
int PsiGrid::i2k( int* i )
{
  int k;

  if (i[1]==0) { // south
      k = i[0];
  } else if (i[0]==Nx) { // east
      k = i[1] + Nx;
  } else if (i[1]==Nx) { // north
      k = 3*Nx - i[0];
  } else if (i[0]==0) { // west
      k = 4*Nx - i[1];
  } else if ( i[0]>0 && i[0]<Nx && i[1]>0 && i[1]<Nx ) { // interior
      k = 4*Nx
        + (Nx-1)*(i[1]-1)
        + i[0]-1;
  } else {
      std::cerr << "Illegal 2D index i=(" << i[0] << "," << i[1] << "). Abort."
                << std::endl;
      exit(EXIT_FAILURE);
  }

  return k;
}
// =============================================================================