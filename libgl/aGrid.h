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
      double Edgelength;
      double h;
      double H0;

      // use the type definition in the CPP file as well (indices!)
      typedef boost::multi_array<double, 2> array_type;
      array_type Ax,
                 Ay;

      float* Ax_i2x( int* );
      float* Ay_i2x( int* );
};