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
    GridVirtual( double scaling = 0.0,
                 double h = 0.0,
                 double gridDomainArea = 0.0,
                 int    numGridPoints = 0,
                 int    numBoundaryPoints = 0);

    virtual
    ~GridVirtual();

    double
    getScaling() const; //!< Returns the scaling factor \f$alpha\f$ of the domain.

    void
    setScaling( const double alpha );

    int
    getNumGridPoints() const; //!< Returns the number of grid points.

    //! Returns the number of grid points on the boundary.
    int
    getNumBoundaryPoints() const;

    double
    getH() const; //!< Returns mesh size \f$h\f$.

    //! Returns the measure of the discretized domain.
    double
    getGridDomainArea() const;

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

    virtual nodeType
    getBoundaryNodeType( int l ) const = 0;

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
    boundaryIndex2globalIndex( int l ) const = 0;

    virtual void
    writeWithGrid( const Tpetra::MultiVector<double,int> & x,
                   const Teuchos::ParameterList &params,
                   const std::string &filePath) const = 0;

  protected:
    double scaling_; //! scaling factor
    double h_;
    double gridDomainArea_;
    unsigned int numGridPoints_;
    unsigned int numBoundaryPoints_;

  private:
  };

  void
  readWithGrid( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
                const std::string                             & filePath,
                Teuchos::RCP<Tpetra::MultiVector<double> >    & x,
                Teuchos::RCP<GridVirtual>                     & grid,
                Teuchos::ParameterList                        & params
              );

#endif /* GRIDVIRTUAL_H_ */
