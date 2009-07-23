#include "aGrid.h"

// include Boost's multidimensional arrays to be able to handle Ax, Ay in a more sane way
#include "boost/multi_array.hpp"

// #include <iostream>

// =============================================================================
// Class constructor
AGrid::AGrid( int nx,
              double edgelength,
              double h0 ):
  Nx(nx),
  Edgelength(edgelength),
  h( edgelength/nx ),
  H0(h0),
  Ax(boost::extents[nx][nx+1]),
  Ay(boost::extents[nx+1][nx])
{
  typedef array_type::index index;

  /*! Initialize the Ax with values
   *  \f[
   *      A_x = - \frac{H_0}{2} y + C.
   *  \f]
   */
  for ( index i=0; i!=nx; ++i )
      for ( index j=0; j!=nx+1; ++j )
          Ax[i][j] = - 0.5 *H0 *j*h
                     + 0.25*H0 *edgelength; //  to level the thing, but not actually necessary

  /*! Initialize the Ay with values
   *  \f[
   *      A_y = \frac{H_0}{2} x + C.
   *  \f]
   */
  for ( index i=0; i!=nx+1; ++i )
      for ( index j=0; j!=nx; ++j )
          Ay[i][j] =   0.5 *H0 *i*h
                     - 0.25*H0 *edgelength; //  to level the thing, but not actually necessary

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
double AGrid::getAxLeft( int *i )
{
  return  Ax[ i[0]-1 ][ i[1] ];
}
// =============================================================================


// =============================================================================
double AGrid::getAxRight( int* i )
{
  return  Ax[ i[0] ][ i[1] ]; // indeed not "+1"; staggered grids!
}
// =============================================================================


// =============================================================================
double AGrid::getAyBelow( int* i )
{
  return  Ay[ i[0] ][ i[1]-1 ];
}
// =============================================================================


// =============================================================================
double AGrid::getAyAbove( int* i )
{
  return  Ay[ i[0] ][ i[1] ]; // indeed not "+1"; staggered grids!
}
// =============================================================================