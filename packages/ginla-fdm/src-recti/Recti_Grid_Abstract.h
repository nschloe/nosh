/*
 * GridVirtual.h
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef RECTI_GRID_ABSTRACT_H_
#define RECTI_GRID_ABSTRACT_H_

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Epetra_MultiVector.h>

#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

#include <LOCA_Parameter_Vector.H>

typedef Tpetra::MultiVector<double              ,Thyra::Ordinal> DoubleMultiVector;

typedef Tpetra::Vector     <std::complex<double>,Thyra::Ordinal> ComplexVector;
typedef Tpetra::MultiVector<std::complex<double>,Thyra::Ordinal> ComplexMultiVector;

typedef Teuchos::Tuple<int,2>           IntTuple;
typedef Teuchos::Tuple<unsigned int,2>  UIntTuple;
typedef Teuchos::Tuple<double,3>        Point;

namespace Recti
{
  namespace Grid
  {

class Abstract
{
public:

    /*!
      Constructor for the grid base class, setting all the essential data at once.
      @param[in] h               Grid distance.
      @param[in] gridDomainArea  Total area occupied by the grid.
      @param[in] numGridPoints   Total number of grid points.
      @param[in] scaling         Scaling of the grid with the given \c h. Defaults to \c 1.0.
    */
    Abstract ( const Point  & h,
               const double         gridDomainArea,
               const unsigned int   numGridPoints,
               const double         scaling=1.0
             );

    //! Empty constructor.
    Abstract();

    virtual
    ~Abstract();

    //! Returns the parameters connected with the grid, e.g, the scaling factor.
    virtual Teuchos::RCP<LOCA::ParameterVector>
    getParameters() const;

    /*!
      Resets the scaling factor of the domain and updates all the dependent values (e.g., \c gridDomainArea_).
      @param[in] alpha New scaling value.
    */
    virtual void
    updateScaling ( const double alpha );

    virtual void
    updateScaling ( const LOCA::ParameterVector & p );

    virtual unsigned int
    getNumGridPoints() const ; //!< Returns the number of grid points.

    virtual Point
    getH() const; //!< Returns mesh sizes \f$h\f$.

    //! Returns the measure of the discretized domain.
    virtual double
    getGridDomainArea() const;

public:
    /*! Indicates whether a node sits in a corner of the domain, on an edge,
     or strictly inside it. */
    enum nodeType
    {
        INTERIOR,  // 0
        BOUNDARY_BOTTOMLEFTCONVEX, // 1
        BOUNDARY_BOTTOMLEFTCONCAVE, // 2
        BOUNDARY_BOTTOMRIGHTCONVEX, // 3
        BOUNDARY_BOTTOMRIGHTCONCAVE, // 4
        BOUNDARY_TOPLEFTCONVEX, // 5
        BOUNDARY_TOPLEFTCONCAVE,  // 6
        BOUNDARY_TOPRIGHTCONVEX, // 7
        BOUNDARY_TOPRIGHTCONCAVE,  // 8
        BOUNDARY_BOTTOM, // 9
        BOUNDARY_RIGHT,  // 10
        BOUNDARY_TOP, // 11
        BOUNDARY_LEFT // 12
    };

public:

    //! @param   k Global node index
    //! @return    Node type
    virtual nodeType
    getNodeType ( unsigned int k ) const = 0;

    //! @param   l Boundary node index
    //! @return    (Boundary) node type
    virtual nodeType
    getBoundaryNodeType ( unsigned int l ) const = 0;

    /*! For a given node number \c k, returns the the area of the surrounding cell. */
    virtual double
    cellArea ( unsigned int k ) const = 0;

    virtual Teuchos::RCP<Point>
    getX ( unsigned int k ) const = 0; //!< Returns the value of \f$x\f$.

    virtual Teuchos::RCP<Point>
    getXLeft ( unsigned int k ) const = 0; //!< Returns the value of \f$x\f$ left of point i.

    virtual Teuchos::RCP<Point>
    getXRight ( unsigned int k ) const = 0; //!< Returns the value of \f$x\f$ right of point i.

    virtual Teuchos::RCP<Point>
    getXBelow ( unsigned int k ) const = 0 ; //!< Returns the value of \f$x\f$ below point i.

    virtual Teuchos::RCP<Point>
    getXAbove ( unsigned int k ) const = 0; //!< Returns the value of \f$x\f$ above point i.

    virtual unsigned int
    getKLeft ( unsigned int k ) const = 0; //!< Returns the running index \c k of the node left of \c i.

    virtual unsigned int
    getKRight ( unsigned int k ) const = 0; //!< Returns the running index \c k of the node right of \c i.

    virtual unsigned int
    getKBelow ( unsigned int k ) const = 0; //!< Returns the running index \c k of the node below \c i.

    virtual unsigned int
    getKAbove ( unsigned int k ) const = 0; //!< Returns the running index \c k of the node above \c i.

    //! Returns the global index \ck of the \cl-th boundary node. Subsequent nodes \c l, \cl+1
    //! sit next to each other.
    virtual unsigned int
    boundaryIndex2globalIndex ( unsigned int l ) const = 0;

    virtual void
    writeWithGrid ( const Epetra_MultiVector     & x,
                    const Teuchos::ParameterList & params,
                    const std::string            & filePath
                  ) const = 0;

    virtual void
    writeWithGrid ( const DoubleMultiVector      & x,
                    const Teuchos::ParameterList & params,
                    const std::string            & filePath
                  ) const = 0;

    virtual void
    writeWithGrid ( const ComplexMultiVector     & x,
                    const Teuchos::ParameterList & params,
                    const std::string            & filePath
                  ) const = 0;

protected:

    Point h_;
    double scaling_; //! scaling factor
    double gridDomainArea_;
    unsigned int numGridPoints_;
    
private:
};

  } // namespace Grid
} // namespace Recti

#endif /* RECTI_GRID_ABSTRACT_H_ */
