/*
 * Grid.h
 *
 *  Created on: Nov 5, 2009
 *      Author: nico
 */

#ifndef GRID_H_
#define GRID_H_

#include <Teuchos_Array.hpp>

class Grid
{
public:

  //! Default constructor.
  Grid(int nx, double edgelength);

  virtual
  ~Grid();

  int
  getNx() const; //!< Returns \f$N_x\f$.

  //! Returns the measure of the discretized domain.
  double
  getGridDomainArea() const;

  double
  getEdgeLength() const; //!< Returns the edge length \f$a\f$ of the square.

  int
  getNumGridPoints() const; //!< Returns the number of grid points.

  double
  getH() const; //!< Returns mesh size \f$h\f$.

  //      int* k2i( int  ); //!< Converts a running index k to a grid index i

  /*! Returns the permutation vector that mats the internal ordering to
   the generic lexicographic ordering of the square grid. */
  void
  lexicographic2grid(std::vector<int> *p) const;

  /*! Indicates whether a node sits in a corner of the domain, on an edge,
   or strictly inside it. */
  enum nodeType
  {
    CORNER, EDGE, INTERIOR
  };

  /*! For a given node number k, indicates whether the node sits in a corner
   of the domain, on an edge, or strictly inside it. */
  nodeType
  k2nodeType(int k) const;

  Teuchos::RCP<Teuchos::Array<double> >
  getXLeft(int k) const; //!< Returns the value of \f$x\f$ left of point i.

  Teuchos::RCP<Teuchos::Array<double> >
  getXRight(int k) const; //!< Returns the value of \f$x\f$ right of point i.

  Teuchos::RCP<Teuchos::Array<double> >
  getXBelow(int k) const; //!< Returns the value of \f$x\f$ below point i.

  Teuchos::RCP<Teuchos::Array<double> >
  getXAbove(int k) const; //!< Returns the value of \f$x\f$ above point i.


  int
  getKLeft(int k) const; //!< Returns the running index \c k of the node left of \c i.

  int
  getKRight(int k) const; //!< Returns the running index \c k of the node right of \c i.

  int
  getKBelow(int k) const; //!< Returns the running index \c k of the node below \c i.

  int
  getKAbove(int k) const; //!< Returns the running index \c k of the node above \c i.

  // TODO: move this to private
  int
  i2k( Teuchos::RCP<Teuchos::Array<int> > & i ) const; //!< Converts a grid index i to a running index k

protected:
private:
  int nx_; //!< Number of grid pieces in both x- and y-direction
  double edgeLength_;
  double h_;

  Teuchos::RCP<Teuchos::Array<double> >
  getX(Teuchos::Array<int> i) const;

  Teuchos::RCP<Teuchos::Array<int> >
  k2i(int k) const;

};

#endif /* GRID_H_ */
