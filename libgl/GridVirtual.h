/*
 * GridVirtual.h
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schlšmer
 */

#ifndef GRIDVIRTUAL_H_
#define GRIDVIRTUAL_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Tpetra_MultiVector.hpp>

class GridVirtual
{
public:
    GridVirtual();

    virtual
    ~GridVirtual();

    //! Returns the measure of the discretized domain.
    virtual double
    getGridDomainArea() const = 0;

    double
    getScaling() const; //!< Returns the scaling factor \f$alpha\f$ of the domain.

    void
    setScaling( const double scaling );

    int
    getNumGridPoints() const; //!< Returns the number of grid points.

    //! Returns the number of grid points on the boundary.
    int
    getNumBoundaryPoints() const;

    double
    getH() const; //!< Returns mesh size \f$h\f$.

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

    nodeType
    boundaryNodeType( int l );

    /*! For a given node number k, returns the the area of the surrounding cell. */
    virtual double
    cellArea(int k) const = 0;

    virtual Teuchos::RCP<Teuchos::Array<double> >
    getXLeft(int k) const = 0; //!< Returns the value of \f$x\f$ left of point i.

    virtual Teuchos::RCP<Teuchos::Array<double> >
    getXRight(int k) const = 0; //!< Returns the value of \f$x\f$ right of point i.

    virtual Teuchos::RCP<Teuchos::Array<double> >
    getXBelow(int k) const = 0 ; //!< Returns the value of \f$x\f$ below point i.

    virtual Teuchos::RCP<Teuchos::Array<double> >
    getXAbove(int k) const = 0; //!< Returns the value of \f$x\f$ above point i.


    virtual int
    getKLeft(int k) const = 0; //!< Returns the running index \c k of the node left of \c i.

    virtual int
    getKRight(int k) const = 0; //!< Returns the running index \c k of the node right of \c i.

    virtual int
    getKBelow(int k) const = 0; //!< Returns the running index \c k of the node below \c i.

    virtual int
    getKAbove(int k) const = 0; //!< Returns the running index \c k of the node above \c i.

    //! Returns the global index \ck of the \cl-th boundary node. Subsequent nodes \c l, \cl+1
    //! sit next to each other.
    virtual int
    boundaryIndex2globalIndex( int l ) = 0;

    // TODO: move this to private
    int
    i2k( Teuchos::RCP<Teuchos::Array<int> > & i ) const; //!< Converts a grid index i to a running index k

    virtual void
    writeWithGrid( const Tpetra::MultiVector<double,int> & x,
                   const Teuchos::ParameterList &params,
                   const std::string &filePath) const = 0;

  protected:
  private:
    double alpha_; //! scaling factor
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
                Teuchos::RCP<GridVirtual>                     & grid,
                Teuchos::ParameterList                        & params
              );

#endif /* GRIDVIRTUAL_H_ */
