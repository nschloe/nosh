/********************************************//**
 * Grid for \f$\psi\f$ with corresponding functions that .
 ***********************************************/
#ifndef STAGGEREDGRID_H
#define STAGGEREDGRID_H

#include <boost/multi_array.hpp>

#include <Teuchos_Array.hpp>

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

      int getNx() const; //!< Returns \f$N_x\f$.

      double getEdgeLength() const; //!< Returns the edge length \f$a\f$ of the square.

      int getNumComplexUnknowns() const; //!< Returns the number of grid points for \f$\psi\f$.

      double getH() const; //!< Returns mesh size \f$h\f$.

//      int* k2i( int  ); //!< Converts a running index k to a grid index i
      int  i2k( Teuchos::Array<int> ) const; //!< Converts a grid index i to a running index k

      /*! Returns the permutation vector that mats the internal ordering to
          the generic lexicographic ordering of the sqare grid. */
      void lexicographic2grid( std::vector<int> *p ) const;

      /*! Indicates whether a node sits in a corner of the domain, on an edge,
          or strictly inside it. */
      enum nodeType { CORNER, EDGE, INTERIOR };

      /*! For a given node number k, indicates whether the node sits in a corner
          of the domain, on an edge, or strictly inside it. */
      nodeType k2nodeType( int k ) const;

      //! Set the magnetic field parameter \f$H_0\f$0.
      void setH0( double h0 );

      double getAxLeft ( Teuchos::Array<int> i ) const;  //!< Returns the value of \f$A_x\f$ left of point i.
      double getAxRight( Teuchos::Array<int> i ) const;  //!< Returns the value of \f$A_x\f$ right of point i.
      double getAyBelow( Teuchos::Array<int> i ) const;  //!< Returns the value of \f$A_y\f$ below point i.
      double getAyAbove( Teuchos::Array<int> i ) const;  //!< Returns the value of \f$A_y\f$ above point i.

      int getKLeft ( Teuchos::Array<int> i ) const; //!< Returns the running index \c k of the node left of \c i.
      int getKRight( Teuchos::Array<int> i ) const; //!< Returns the running index \c k of the node right of \c i.
      int getKBelow( Teuchos::Array<int> i ) const; //!< Returns the running index \c k of the node below \c i.
      int getKAbove( Teuchos::Array<int> i ) const; //!< Returns the running index \c k of the node above \c i.

  private:
      int    nx_; //!< Number of grid pieces in both x- and y-direction
      double edgeLength_;
      double h0_;
      double h_;

      // use the type definition in the CPP file as well (indices!)
      typedef boost::multi_array<double,2> array_type;
      array_type Ax_;
      array_type Ay_;

      //! (Re)sets the values of \f$A\f$.
      void computeA();
};
#endif // STAGGEREDGRID_H