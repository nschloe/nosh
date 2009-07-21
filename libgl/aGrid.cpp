#include "aGrid.h"

// include Boost's multidimensional arrays to be able to handle Ax, Ay in a more sane way
#include "boost/multi_array.hpp"

#include <iostream>

// =============================================================================
// Class constructor
AGrid::AGrid( int nx,
              double edgelength,
              double h0 ):
Nx(nx),
Edgelength(edgelength),
h(0.0),
H0(h0),
Ax(boost::extents[nx][nx+1]),
Ay(boost::extents[nx+1][nx])
{
  // initialize the step size
  h = Edgelength / Nx;

  typedef array_type::index index;

  // initialize the Ax with values
  for ( index i=0; i!=nx; ++i )
      for ( index j=0; j!=nx+1; ++j )
          Ax[i][j] = - 0.5 *H0 *j*h
                     + 0.25*H0 *edgelength; //  to level the thing, but not actually necessary

  // initialize the Ay with values
  for ( index i=0; i!=nx+1; ++i )
      for ( index j=0; j!=nx; ++j )
          Ay[i][j] =   0.5 *H0 *i*h
                     + 0.25*H0 *edgelength; //  to level the thing, but not actually necessary

//   // ---------------------------------------------------------------------------
//   // for debugging purposes:
//   std::cout << "Nx=" << Nx << std::endl;
//   std::cout << "Edgelength=" << Edgelength << std::endl;
//   std::cout << "h=" << h << std::endl;
//   std::cout << "H0=" << H0 << std::endl;
//   for ( index i=0; i!=nx; ++i )
//       for ( index j=0; j!=nx+1; ++j )
//           std::cout << "Ax[" << i << "][" << j << "] = " << Ax[i][j]
//                     << std::endl;
//   // ---------------------------------------------------------------------------

}
// =============================================================================


// =============================================================================
// Destructor
AGrid::~AGrid()
{
}
// =============================================================================


// =============================================================================
// converts an integer index to the geometric positio
float* AGrid::Ax_i2x( int* i )
{
  static float x[2];

  x[0] = i[0]*h + h/2;
  x[1] = i[1]*h;

  return x;
}
// =============================================================================



// =============================================================================
// converts an integer index to the geometric positio
float* AGrid::Ay_i2x( int* i )
{
  static float x[2];

  x[0] = i[0]*h;
  x[1] = i[1]*h + h/2;

  return x;
}
// =============================================================================



// =============================================================================
float AGrid::getAxLeft( int *i )
{
  return  Ax[ i[0]-1 ][ i[1] ];
}
// =============================================================================


// =============================================================================
float AGrid::getAxRight( int* i )
{
  return  Ax[ i[0] ][ i[1] ]; // indeed not "+1"; staggered grids!
}
// =============================================================================


// =============================================================================
float AGrid::getAyBelow( int* i )
{
  return  Ay[ i[0] ][ i[1]-1 ];
}
// =============================================================================


// =============================================================================
float AGrid::getAyAbove( int* i )
{
  return  Ay[ i[0] ][ i[1] ]; // indeed not "+1"; staggered grids!
}
// =============================================================================