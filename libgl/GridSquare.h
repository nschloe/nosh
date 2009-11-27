/*
 * Grid.h
 *
 *  Created on: Nov 5, 2009
 *      Author: Nico Schlšmer
 */

#ifndef GRIDSQUARE_H_
#define GRIDSQUARE_H_

#include "GridVirtual.h"

#include <Teuchos_Array.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Comm.hpp>

typedef Tpetra::Vector<double>      DoubleVector;
typedef Tpetra::MultiVector<double> DoubleMultiVector;

class GridSquare: public GridVirtual
{

public:

  //! Default constructor.
  GridSquare( int nx = 0,
              double scaling = 0);

  virtual
  ~GridSquare();

  virtual nodeType
  getBoundaryNodeType( int l ) const;

  /*! For a given node number k, returns the the area of the surrounding cell. */
  virtual double
  cellArea(int k) const;

  virtual Teuchos::RCP<Teuchos::Array<double> >
  getXLeft(int k) const; //!< Returns the value of \f$x\f$ left of point i.

  virtual Teuchos::RCP<Teuchos::Array<double> >
  getXRight(int k) const; //!< Returns the value of \f$x\f$ right of point i.

  virtual Teuchos::RCP<Teuchos::Array<double> >
  getXBelow(int k) const; //!< Returns the value of \f$x\f$ below point i.

  virtual Teuchos::RCP<Teuchos::Array<double> >
  getXAbove(int k) const; //!< Returns the value of \f$x\f$ above point i.


  virtual int
  getKLeft(int k) const; //!< Returns the running index \c k of the node left of \c i.

  virtual int
  getKRight(int k) const; //!< Returns the running index \c k of the node right of \c i.

  virtual int
  getKBelow(int k) const; //!< Returns the running index \c k of the node below \c i.

  virtual int
  getKAbove(int k) const; //!< Returns the running index \c k of the node above \c i.

  //! Returns the global index \ck of the \cl-th boundary node. Subsequent nodes \c l, \cl+1
  //! sit next to each other.
  virtual int
  boundaryIndex2globalIndex( int l ) const;

  void
  writeWithGrid( const Tpetra::MultiVector<double,int> & x,
                 const Teuchos::ParameterList &params,
                 const std::string &filePath) const;

  void
  read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
        const std::string                             & filePath,
        Teuchos::RCP<DoubleMultiVector>               & x,
        Teuchos::ParameterList                        & params
      );

protected:
private:
  int nx_; //!< Number of grid pieces in both x- and y-direction

  Teuchos::RCP<DoubleMultiVector>
  permuteGrid2Lexicographic( const DoubleMultiVector & x ) const;

  Teuchos::RCP<DoubleMultiVector>
  permuteLexicographic2Grid( const DoubleMultiVector & xLexicographic
                           ) const;

  //! Defines a subsequent order of boundary nodes by associating a running index \c l with
  //! \f$i\f$-coordinates on the grid.
  Teuchos::RCP<Teuchos::Array<int> >
  boundaryPosition ( int l ) const;

  Teuchos::RCP<Teuchos::Array<double> >
  getX(Teuchos::Array<int> i) const;

  Teuchos::RCP<Teuchos::Array<int> >
  k2i(int k) const;

  int
  i2k( Teuchos::RCP<Teuchos::Array<int> > & i ) const; //!< Converts a grid index i to a running index k

};

#endif /* GRIDSQUARE_H_ */
