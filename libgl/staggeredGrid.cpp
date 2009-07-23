#include <iostream>
#include "staggeredGrid.h"

// =============================================================================
// Class constructor
StaggeredGrid::StaggeredGrid( int    nx,
                              double edgelength,
                              double h0          ):
  Nx(nx),
  Edgelength(edgelength),
  H0(h0),
  h( edgelength/nx ),
  Ax(boost::extents[nx][nx+1]),
  Ay(boost::extents[nx+1][nx])
{
  // setup A
  computeA();
}
// =============================================================================


// =============================================================================
// Destructor
StaggeredGrid::~StaggeredGrid()
{
}
// =============================================================================


// =============================================================================
int StaggeredGrid::getNx()
{
    return Nx;
}
// =============================================================================


// =============================================================================
int StaggeredGrid::getNumComplexUnknowns()
{
    // the number of grid points in psi
    return (Nx+1)*(Nx+1);
}
// =============================================================================


// =============================================================================
double StaggeredGrid::getH()
{
    return h;
}
// =============================================================================


// =============================================================================
void StaggeredGrid::setH0( double h0 )
{
  H0 = h0;
  // rebuild the magnetic vector potential values
  computeA();
}
// =============================================================================




// =============================================================================
void StaggeredGrid::computeA()
{
  typedef array_type::index index;

  /*! Initialize the Ax with values
   *  \f[
   *      A_x = - \frac{H_0}{2} y + C.
   *  \f]
   */
  for ( index i=0; i!=Nx; ++i )
      for ( index j=0; j!=Nx+1; ++j )
          Ax[i][j] = - 0.5 *H0 *j*h
                     + 0.25*H0 *Edgelength; //  to level the thing, but not actually necessary

  /*! Initialize the Ay with values
   *  \f[
   *      A_y = \frac{H_0}{2} x + C.
   *  \f]
   */
  for ( index i=0; i!=Nx+1; ++i )
      for ( index j=0; j!=Nx; ++j )
          Ay[i][j] =   0.5 *H0 *i*h
                     - 0.25*H0 *Edgelength; //  to level the thing, but not actually necessary

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
double StaggeredGrid::getAxLeft( int *i )
{
  return  Ax[ i[0]-1 ][ i[1] ];
}
// =============================================================================


// =============================================================================
double StaggeredGrid::getAxRight( int* i )
{
  return  Ax[ i[0] ][ i[1] ]; // indeed not "+1"; staggered grids!
}
// =============================================================================


// =============================================================================
double StaggeredGrid::getAyBelow( int* i )
{
  return  Ay[ i[0] ][ i[1]-1 ];
}
// =============================================================================


// =============================================================================
double StaggeredGrid::getAyAbove( int* i )
{
  return  Ay[ i[0] ][ i[1] ]; // indeed not "+1"; staggered grids!
}
// =============================================================================


// // =============================================================================
// // maps a running index k to a 2D index i
// int* StaggeredGrid::k2i( int k )
// {
//   static int i[2];
//
//   if (k<Nx) { // lower shore
//     i[0] = k-1;
//     i[1] = 0;
//   }
//   else if (k<2*Nx) { // right shore
//     i[0] = Nx;
//     i[1] = k-Nx;
//   }
//   else if (k<3*Nx) { // upper shore
//     i[0] = 3*Nx-k;
//     i[1] = Nx;
//   }
//   else if (k<4*Nx) { // left shore
//     i[0] = 0;
//     i[1] = 4*Nx-k;
//   }
//   else { // on the interior
//     int numBoundaryNodes = 4*Nx;
//     i[0] = (k-numBoundaryNodes)%(Nx-1) + 1;
//     i[1] = (k-numBoundaryNodes)/(Nx-1) + 1;
//   }
//
//   return i;
// }
// // =============================================================================


// =============================================================================
StaggeredGrid::nodeType StaggeredGrid::k2nodeType( int k )
{
  if (k==0 || k==Nx || k==2*Nx || k==3*Nx )
      return StaggeredGrid::CORNER;
  else if (k<4*Nx)
      return StaggeredGrid::EDGE;
  else
      return StaggeredGrid::INTERIOR;
}
// =============================================================================


// =============================================================================
// maps a 2D index i to a running index k
int StaggeredGrid::i2k( int* i )
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
      std::cerr << "ERROR: PsiGrid::i2k - "
                << "    Illegal 2D index i=(" << i[0] << "," << i[1] << ")."
                << " Abort."
                << std::endl;
      exit(EXIT_FAILURE);
  }

  return k;
}
// =============================================================================