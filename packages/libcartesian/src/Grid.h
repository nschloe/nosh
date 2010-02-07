/*
 * Grid.h
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef GRID_H_
#define GRID_H_

#include "DomainVirtual.h"
#include "GridVirtual.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Epetra_MultiVector.h>

#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal


class Grid:
            virtual public GridVirtual
{
public:

    //! Default constructor.
    Grid ( const Teuchos::RCP<const DomainVirtual> & domain,
           const DoubleTuple                         h,
           const double                              scaling = 1.0
         );

    //! Empty constructor.
    Grid();

    virtual
    ~Grid();

    virtual double
    getScaling() const; //!< Returns the scaling factor \f$alpha\f$ of the domain.

    virtual void
    setScaling ( const double alpha );

    virtual unsigned int
    getNumGridPoints() const ; //!< Returns the number of grid points.

    //! Returns the number of grid points on the boundary.
    virtual unsigned int
    getNumBoundaryPoints() const;

    virtual DoubleTuple
    getH() const; //!< Returns mesh sizes \f$h\f$.

    //! Returns the measure of the discretized domain.
    virtual double
    getGridDomainArea() const;

public:

    //! @param   k Global node index
    //! @return    Node type
    virtual nodeType
    getNodeType ( unsigned int k ) const;

    //! @param   l Boundary node index
    //! @return    (Boundary) node type
    virtual nodeType
    getBoundaryNodeType ( unsigned int l ) const;

    /*! For a given node number \c k, returns the the area of the surrounding cell. */
    virtual double
    cellArea ( unsigned int k ) const;

    virtual Teuchos::RCP<DoubleTuple>
    getXLeft ( unsigned int k ) const; //!< Returns the value of \f$x\f$ left of point i.

    Teuchos::RCP<DoubleTuple>
    getXLeft ( const IntTuple & i ) const;

    virtual Teuchos::RCP<DoubleTuple>
    getXRight ( unsigned int k ) const; //!< Returns the value of \f$x\f$ right of point i.

    Teuchos::RCP<DoubleTuple>
    getXRight ( const IntTuple & i ) const;

    virtual Teuchos::RCP<DoubleTuple>
    getXBelow ( unsigned int k ) const ; //!< Returns the value of \f$x\f$ below point i.

    Teuchos::RCP<DoubleTuple>
    getXBelow ( const IntTuple & i ) const;

    virtual Teuchos::RCP<DoubleTuple>
    getXAbove ( unsigned int k ) const; //!< Returns the value of \f$x\f$ above point i.

    Teuchos::RCP<DoubleTuple>
    getXAbove ( const IntTuple & i ) const;

    virtual unsigned int
    getKLeft ( unsigned int k ) const; //!< Returns the running index \c k of the node left of \c i.

    virtual unsigned int
    getKRight ( unsigned int k ) const; //!< Returns the running index \c k of the node right of \c i.

    virtual unsigned int
    getKBelow ( unsigned int k ) const; //!< Returns the running index \c k of the node below \c i.

    virtual unsigned int
    getKAbove ( unsigned int k ) const; //!< Returns the running index \c k of the node above \c i.

    //! Returns the global index \ck of the \cl-th boundary node. Subsequent nodes \c l, \cl+1
    //! sit next to each other.
    virtual unsigned int
    boundaryIndex2globalIndex ( unsigned int l ) const;

    virtual void
    writeWithGrid ( const Epetra_MultiVector     & x,
                    const Teuchos::ParameterList & params,
                    const std::string            & filePath
                  ) const;

    virtual void
    writeWithGrid ( const DoubleMultiVector      & x,
                    const Teuchos::ParameterList & params,
                    const std::string            & filePath
                  ) const;

    virtual void
    writeWithGrid ( const ComplexMultiVector     & x,
                    const Teuchos::ParameterList & params,
                    const std::string            & filePath
                  ) const;

    void
    writeWithBoundingBoxGrid ( const Epetra_MultiVector     & x,
                               const Teuchos::ParameterList & params,
                               const std::string            & filePath
                             ) const;

    virtual void
    read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
           const std::string                             & filePath,
           Teuchos::RCP<DoubleMultiVector>               & x,
           Teuchos::ParameterList                        & params
         );

    virtual void
    read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
           const std::string                             & filePath,
           Teuchos::RCP<ComplexMultiVector>              & x,
           Teuchos::ParameterList                        & params
         );


private:

    const Teuchos::RCP<const DomainVirtual> domain_;
    DoubleTuple h_;
    UIntTuple Nx_;
    Teuchos::Array<int> kBB_;
    double scaling_; //! scaling factor
    double gridDomainArea_;
    unsigned int numGridPoints_;
    unsigned int numBoundaryPoints_;
    Teuchos::Array<IntTuple> nodes_;
    Teuchos::Array<nodeType> nodeTypes_;
    Teuchos::Tuple<double,2> origin_;

private:
    //! Indicates directions when traversing the boundary of a
    //! domain.
    enum direction
    {
        LEFT,  // 0
        RIGHT, // 1 
        UP,    // 2
        DOWN   // 3
    };

private:
    Grid::nodeType
    getNodeType ( const direction prevDir,
                  const direction newDir
                ) const;

    Teuchos::RCP<DoubleTuple>
    getX ( const IntTuple & i ) const;

    bool
    boundaryStepper ( Teuchos::Array<IntTuple>  & boundaryNodes,
                      Teuchos::Array<direction> & directions
                    ) const;

    bool
    equal ( const IntTuple & a,
            const IntTuple & b
          ) const;

    IntTuple
    findFirstBoundaryNode () const;

    Grid::direction
    findNextDirection ( const IntTuple  & currentBoundaryNode,
                        const direction   prevDir
                      ) const;

    Teuchos::Tuple<Grid::direction,3>
    getNextDirections ( const direction dir ) const;

    IntTuple
    step ( const IntTuple  & node,
           const direction   dir
         ) const;

    unsigned int
    i2kBoundingBox ( const IntTuple & i ) const;

    Teuchos::RCP<IntTuple>
    k2iBoundingBox ( const unsigned int k ) const;

    unsigned int
    kBoundingBox2kDomain ( const unsigned int k ) const;

    void
    pruneInitialTentacle ( Teuchos::Array<IntTuple>  & nodes,
                           Teuchos::Array<direction> & directions
                         ) const;
};

#endif /* GRID_H_ */
