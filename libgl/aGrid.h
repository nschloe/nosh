// include Boost's multidimensional arrays to be able to handle Ax, Ay in a more sane way
#include "boost/multi_array.hpp"

class AGrid
{
  public:

     AGrid( int nx,
            double edgelength,
            double h0 );

     ~AGrid();

      float getAxLeft ( int* );
      float getAxRight( int* );
      float getAyBelow( int* );
      float getAyAbove( int* );

//   protected:
//      // Attributes visible to descendents
  private:
      int Nx;
      int Edgelength;
      double h;
      int d;
      double H0;
      boost::multi_array<double, 2> Ax;
      boost::multi_array<double, 2> Ay;

      float* Ax_i2x( int* );
      float* Ay_i2x( int* );
};