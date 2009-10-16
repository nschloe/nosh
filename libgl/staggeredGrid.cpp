#include "staggeredGrid.h"
#include "glException.h"

#include <iostream>

#include <EpetraExt_Utils.h>

// =============================================================================
// Class constructor
StaggeredGrid::StaggeredGrid( int    nx,
                              double edgelength,
                              double h0          ):
  nx_(nx),
  edgeLength_(edgelength),
  h0_(h0),
  h_( edgelength/nx ),
  Ax_(boost::extents[nx][nx+1]),
  Ay_(boost::extents[nx+1][nx])
{
  // setup A
  computeA();
}
// =============================================================================
// copy constructor
StaggeredGrid::StaggeredGrid(const StaggeredGrid& sGrid) {
    nx_         = sGrid.nx_;
    edgeLength_ = sGrid.edgeLength_;
    h0_         = sGrid.h0_;
    h_          = sGrid.h_;

    Ax_.resize(boost::extents[nx_][nx_+1]);
    Ax_ = sGrid.Ax_;
    Ay_.resize(boost::extents[nx_+1][nx_]);
    Ay_ = sGrid.Ay_;
}
// =============================================================================
// Destructor
StaggeredGrid::~StaggeredGrid()
{
}
// =============================================================================
int StaggeredGrid::getNx() const
{
    return nx_;
}
// =============================================================================
double StaggeredGrid::getEdgeLength() const
{
    return edgeLength_;
}
// =============================================================================
int StaggeredGrid::getNumComplexUnknowns() const
{
    // the number of grid points in psi
    return (nx_+1)*(nx_+1);
}
// =============================================================================
double StaggeredGrid::getH() const
{
    return h_;
}
// =============================================================================
void StaggeredGrid::setH0( double h0 )
{
  h0_ = h0;
  // rebuild the magnetic vector potential values
  computeA();
}
// =============================================================================
void StaggeredGrid::computeA()
{
  typedef array_type::index index;

  /*! Initialize the Ax_ with values
   *  \f[
   *      A_x = - \frac{H_0}{2} y + C.
   *  \f]
   */
  for ( index i=0; i!=nx_; ++i )
      for ( index j=0; j!=nx_+1; ++j )
          Ax_[i][j] = - 0.5*h0_ *j*h_
                      +0.25*h0_ *edgeLength_; //  to level the thing, but not actually necessary

  /*! Initialize the Ay_ with values
   *  \f[
   *      A_y = \frac{H_0}{2} x + C.
   *  \f]
   */
  for ( index i=0; i!=nx_+1; ++i )
      for ( index j=0; j!=nx_; ++j )
          Ay_[i][j] =   0.5*h0_ *i*h_
                     - 0.25*h0_ *edgeLength_; //  to level the thing, but not actually necessary

//   // ---------------------------------------------------------------------------
//   // for debugging purposes:
//   std::cout << "nx_=" << nx_ << std::endl;
//   std::cout << "edgeLength_=" << edgeLength_ << std::endl;
//   std::cout << "h=" << h << std::endl;
//   std::cout << "h0_=" << h0_ << std::endl;
//   for ( index i=0; i!=nx_; ++i )
//       for ( index j=0; j!=nx_+1; ++j )
//           std::cout << "Ax_[" << i << "][" << j << "] = " << Ax_[i][j]
//                     << std::endl;
//   // ---------------------------------------------------------------------------

}
// =============================================================================
double StaggeredGrid::getAxLeft( Teuchos::Array<int> i ) const
{
  return  Ax_[ i[0]-1 ][ i[1] ];
}
// =============================================================================
double StaggeredGrid::getAxRight( Teuchos::Array<int> i ) const
{
  return  Ax_[ i[0] ][ i[1] ]; // indeed not "+1"; staggered grids!
}
// =============================================================================
double StaggeredGrid::getAyBelow( Teuchos::Array<int> i ) const
{
  return  Ay_[ i[0] ][ i[1]-1 ];
}
// =============================================================================
double StaggeredGrid::getAyAbove( Teuchos::Array<int> i ) const
{
  return  Ay_[ i[0] ][ i[1] ]; // indeed not "+1"; staggered grids!
}
// =============================================================================
int StaggeredGrid::getKLeft( Teuchos::Array<int> i ) const
{
  Teuchos::Array<int> j = Teuchos::tuple( i[0]-1, i[1] );
  return i2k( j );
}
// =============================================================================
int StaggeredGrid::getKRight( Teuchos::Array<int> i ) const
{
  Teuchos::Array<int> j = Teuchos::tuple( i[0]+1, i[1] );
  return i2k( j );
}
// =============================================================================
int StaggeredGrid::getKBelow( Teuchos::Array<int> i ) const
{
  Teuchos::Array<int> j = Teuchos::tuple( i[0], i[1]-1 );
  return i2k( j );
}
// =============================================================================
int StaggeredGrid::getKAbove( Teuchos::Array<int> i ) const
{
  Teuchos::Array<int> j = Teuchos::tuple( i[0], i[1]+1 );
  return i2k( j );
}
// =============================================================================
// // maps a running index k to a 2D index i
// int* StaggeredGrid::k2i( int k )
// {
//   static int i[2];
//
//   if (k<nx_) { // lower shore
//     i[0] = k-1;
//     i[1] = 0;
//   }
//   else if (k<2*nx_) { // right shore
//     i[0] = nx_;
//     i[1] = k-nx_;
//   }
//   else if (k<3*nx_) { // upper shore
//     i[0] = 3*nx_-k;
//     i[1] = nx_;
//   }
//   else if (k<4*nx_) { // left shore
//     i[0] = 0;
//     i[1] = 4*nx_-k;
//   }
//   else { // on the interior
//     int numBoundaryNodes = 4*nx_;
//     i[0] = (k-numBoundaryNodes)%(nx_-1) + 1;
//     i[1] = (k-numBoundaryNodes)/(nx_-1) + 1;
//   }
//
//   return i;
// }
// // =============================================================================
StaggeredGrid::nodeType StaggeredGrid::k2nodeType( int k ) const
{
  if (k==0 || k==nx_ || k==2*nx_ || k==3*nx_ )
      return StaggeredGrid::CORNER;
  else if (k<4*nx_)
      return StaggeredGrid::EDGE;
  else
      return StaggeredGrid::INTERIOR;
}
// =============================================================================
// maps a 2D index i to a running index k
int StaggeredGrid::i2k( Teuchos::Array<int> i ) const
{
  int k;

  if (i[1]==0) { // south
      k = i[0];
  } else if (i[0]==nx_) { // east
      k = i[1] + nx_;
  } else if (i[1]==nx_) { // north
      k = 3*nx_ - i[0];
  } else if (i[0]==0) { // west
      k = 4*nx_ - i[1];
  } else if ( i[0]>0 && i[0]<nx_ && i[1]>0 && i[1]<nx_ ) { // interior
      k = 4*nx_
        + (nx_-1)*(i[1]-1)
        + i[0]-1;
  } else {
      std::string message = "Illegal 2D index   i = " + Teuchos::toString( i );
      throw glException( "StaggeredGrid::i2k",
                         message );
  }

  return k;
}
// =============================================================================
// Returns a vector that defines the reordering from a lexicographic grid to
// the ordering present in this grid
void StaggeredGrid::lexicographic2grid( std::vector<int> *p ) const
{
  // check if for admissible vector size
  unsigned int numUnknowns = (nx_+1)*(nx_+1);
  if ( p->size() != numUnknowns ) {
      std::string message = "Size of the input vector p ("
                          + EpetraExt::toString( int(p->size()) ) + ") "
                          + "does not coincide with with number of unknowns on "
                          + " the grid (" + EpetraExt::toString( (nx_+1)*(nx_+1) )
                          + ").";
      throw glException( "StaggeredGrid::lexicographic2grid",
                         message );
  }

  int k=0;
  Teuchos::Array<int> index(2);
  for (int j=0; j<nx_+1; j++) {
      index[1]  = j;
      for (int i=0; i<nx_+1; i++) {
          index[0] = i;
          (*p)[k++] = i2k(index);
      }
  }

}
// =============================================================================