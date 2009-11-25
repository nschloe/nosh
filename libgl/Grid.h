/*
 * Grid.h
 *
 *  Created on: Nov 5, 2009
 *      Author: Nico Schl�mer
 */

#ifndef GRID_H_
#define GRID_H_

#include <Teuchos_Array.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Comm.hpp>

class Grid
{
public:

  //! Default constructor.
  Grid(int nx, double edgeLength);

  virtual
  ~Grid();

  int
  getNx() const; //!< Returns \f$N_x\f$.

  //! Returns the measure of the discretized domain.
  double
  getGridDomainArea() const;

  double
  getEdgeLength() const; //!< Returns the edge length \f$a\f$ of the square.

  void
  setEdgeLength( const double edgeLength );

  int
  getNumGridPoints() const; //!< Returns the number of grid points.

  //! Returns the number of grid points on the boundary.
  int
  getNumBoundaryPoints() const;

  double
  getH() const; //!< Returns mesh size \f$h\f$.

  /*! Returns the permutation vector that mats the internal ordering to
   the generic lexicographic ordering of the square grid. */
  void
  lexicographic2grid(std::vector<int> *p) const;

  /*! Indicates whether a node sits in a corner of the domain, on an edge,
   or strictly inside it. */
  enum nodeType
  {
    BOTTOMLEFTCONVEX,
    BOTTOMLEFTCONCAVE,
    BOTTOMRIGHTCONVEX,
    BOTTOMRIGHTCONCAVE,
    TOPLEFTCONVEX,
    TOPLEFTCONCAVE,
    TOPRIGHTCONVEX,
    TOPRIGHTCONCAVE,
    BOTTOM,
    RIGHT,
    TOP,
    LEFT,
    INTERIOR
  };

  Grid::nodeType
  boundaryNodeType( int l );

  /*! For a given node number k, returns the the area of the surrounding cell. */
  double
  cellArea(int k) const;

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

  //! Returns the global index \ck of the \cl-th boundary node. Subsequent nodes \c l, \cl+1
  //! sit next to each other.
  int
  boundaryIndex2globalIndex( int l );

  // TODO: move this to private
  int
  i2k( Teuchos::RCP<Teuchos::Array<int> > & i ) const; //!< Converts a grid index i to a running index k

  void
  reorderToLexicographic( Tpetra::Vector<std::complex<double> > & x
                        ) const;

  void
  reorderFromLexicographic( Tpetra::Vector<std::complex<double> > & x
                          ) const;

  void
  writeWithGrid( const Tpetra::MultiVector<double,int> & x,
                 const Teuchos::ParameterList &params,
                 const std::string &filePath) const;

protected:
private:
  int nx_; //!< Number of grid pieces in both x- and y-direction
  double edgeLength_;
  double h_;

  //! Defines a subsequent order of boundary nodes by associating a running index \c l with
  //! \f$i\f$-coordinates on the grid.
  Teuchos::RCP<Teuchos::Array<int> >
  boundaryPosition ( int l );

  Teuchos::RCP<Teuchos::Array<double> >
  getX(Teuchos::Array<int> i) const;

  Teuchos::RCP<Teuchos::Array<int> >
  k2i(int k) const;

};

void
readWithGrid( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
              const std::string                             & filePath,
              Teuchos::RCP<Tpetra::MultiVector<double> >    & x,
              Teuchos::RCP<Grid>                            & grid,
              Teuchos::ParameterList                        & params
            );

#endif /* GRID_H_ */
