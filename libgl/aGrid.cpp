#include "aGrid.h"

// include Boost's multidimensional arrays to be able to handle Ax, Ay in a more sane way
#include "boost/multi_array.hpp"


// =============================================================================
// Class constructor
AGrid::AGrid( int nx,
              double edgelength,
              double h0 ):
Nx(nx),
Edgelength(edgelength),
h(0.0),
d(2),
H0(h0),
Ax(boost::extents[nx][nx+1]),
Ay(boost::extents[nx+1][nx])
{
  h = Edgelength / Nx;

  // initialize the Ax with values
  for ( int i=0; i<nx; i++ )
      for ( int j=0; j<nx+1; i++ )
          Ax[i][j] = -h0/2 * j*h;

  // initialize the Ay with values
  for ( int i=0; i<nx+1; i++ )
      for ( int j=0; j<nx; i++ )
          Ay[i][j] = h0/2 * i*h;

}
// =============================================================================


// =============================================================================
// Destructor
AGrid::~AGrid()
{
// delete the grids
}
// =============================================================================


// =============================================================================
// converts an integer index to the geometric positio
float* AGrid::Ax_i2x( int* i )
{
  static float x[2];

  x[1] = i[1]*h + h/2;
  x[2] = i[2]*h;

  return x;
}
// =============================================================================



// =============================================================================
// converts an integer index to the geometric positio
float* AGrid::Ay_i2x( int* i )
{
  static float x[2];

  x[1] = i[1]*h;
  x[2] = i[2]*h + h/2;

  return x;
}
// =============================================================================



// =============================================================================
float AGrid::getAxLeft( int *i )
{
  return  Ax[ i[1]-1 ][ i[2] ];
}
// =============================================================================


// =============================================================================
float AGrid::getAxRight( int* i )
{
  return  Ax[ i[1] ][ i[2] ];
}
// =============================================================================


// =============================================================================
float AGrid::getAyBelow( int* i )
{
  return  Ay[ i[1] ][ i[2]-1 ];
}
// =============================================================================


// =============================================================================
float AGrid::getAyAbove( int* i )
{
  return  Ay[ i[1] ][ i[2] ];
}
// =============================================================================