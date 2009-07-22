/********************************************//**
 * Grid for \f$A\f$.
 ***********************************************/

// include Boost's multidimensional arrays to be able to handle Ax, Ay in a more sane way
#include "boost/multi_array.hpp"

class AGrid
{
  public:

     //! Class constructor.
     AGrid( int nx,
            double edgelength,
            double h0 );

     //! Class destructor
     ~AGrid();

      double getAxLeft ( int* i );  //!< Returns the value of \f$A_x\f$ left of point i.
      double getAxRight( int* i );  //!< Returns the value of \f$A_x\f$ right of point i.
      double getAyBelow( int* i );  //!< Returns the value of \f$A_y\f$ below point i.
      double getAyAbove( int* i );  //!< Returns the value of \f$A_y\f$ above point i.

  private:
      int Nx;
      double Edgelength;
      double h;
      double H0;

      // use the type definition in the CPP file as well (indices!)
      typedef boost::multi_array<double, 2> array_type;
      array_type Ax,
                 Ay;
};