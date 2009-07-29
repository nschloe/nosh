/********************************************//**
 * Grid for \f$\psi\f$ with corresponding functions that .
 ***********************************************/
#ifndef STAGGEREDGRID_H
#define STAGGEREDGRID_H


#include "boost/multi_array.hpp"

class StaggeredGrid
{
  public:
      //! Default constructor.
      StaggeredGrid( int    nx,
                     double edgelength,
                     double h0 );

      //! Copy constructor.
      StaggeredGrid(const StaggeredGrid& sGrid);

      ~StaggeredGrid();

      int getNx(); //!< Returns \f$N_x\f$.

      int getNumComplexUnknowns(); //!< Returns the number of grid points for \f$\psi\f$.

      double getH(); //!< Returns mesh size \f$h\f$.

//      int* k2i( int  ); //!< Converts a running index k to a grid index i
      int  i2k( int* ); //!< Converts a grid index i to a running index k

      /*! Indicates whether a node sits in a corner of the domain, on an edge,
          or strictly inside it. */
      enum nodeType { CORNER, EDGE, INTERIOR };

      /*! For a given node number k, indicates whether the node sits in a corner
          of the domain, on an edge, or strictly inside it. */
      nodeType k2nodeType( int k );

      //! Set the magnetic field parameter \f$H_0\f$0.
      void setH0( double h0 );

      double getAxLeft ( int* i );  //!< Returns the value of \f$A_x\f$ left of point i.
      double getAxRight( int* i );  //!< Returns the value of \f$A_x\f$ right of point i.
      double getAyBelow( int* i );  //!< Returns the value of \f$A_y\f$ below point i.
      double getAyAbove( int* i );  //!< Returns the value of \f$A_y\f$ above point i.

  private:
      int    Nx; //!< Number of grid pieces in both x- and y-direction
      double Edgelength,
             H0,
             h;

      // use the type definition in the CPP file as well (indices!)
      typedef boost::multi_array<double,2> array_type;
      array_type Ax,
                 Ay;

      //! (Re)sets the values of \f$A\f$.
      void computeA();
};
#endif // STAGGEREDGRID_H