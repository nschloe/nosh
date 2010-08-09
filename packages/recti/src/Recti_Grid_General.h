/*
 * Grid.h
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef GRID_GENERAL_H_
#define GRID_GENERAL_H_

#include "Recti_Domain_Abstract.h"
#include "Recti_Grid_Abstract.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Epetra_MultiVector.h>

#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

namespace Recti
{
  namespace Grid
  {

class General:
            public Abstract
{
public:

    //! Default constructor.
    General ( const Teuchos::RCP<const Recti::Domain::Abstract> & domain,
              const Point                                 & h
            );
         
    //! Constructor for information gotten from a file read.
    General ( const Point         & h,
              const UIntTuple           & numCells,
              const Teuchos::Array<int> & kBB,
              const Teuchos::Array<int> & boundaryNodes,
              const double                scaling,
              const Point         & origin
            );

    //! Empty constructor.
    General();

    virtual
    ~General();
    
    virtual const Teuchos::Array<int> &
    getBoundaryIndices() const;

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

    virtual Teuchos::RCP<Point>
    getX ( unsigned int k ) const; //!< Returns the value of \f$x\f$. */

    virtual Teuchos::RCP<Point>
    getXLeft ( unsigned int k ) const; //!< Returns the value of \f$x\f$ left of point i.

    Teuchos::RCP<Point>
    getXLeft ( const UIntTuple & i ) const;

    virtual Teuchos::RCP<Point>
    getXRight ( unsigned int k ) const; //!< Returns the value of \f$x\f$ right of point i.

    Teuchos::RCP<Point>
    getXRight ( const UIntTuple & i ) const;

    virtual Teuchos::RCP<Point>
    getXBelow ( unsigned int k ) const ; //!< Returns the value of \f$x\f$ below point i.

    Teuchos::RCP<Point>
    getXBelow ( const UIntTuple & i ) const;

    virtual Teuchos::RCP<Point>
    getXAbove ( unsigned int k ) const; //!< Returns the value of \f$x\f$ above point i.

    Teuchos::RCP<Point>
    getXAbove ( const UIntTuple & i ) const;

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

private:

    UIntTuple numCells_;
    Teuchos::Array<int> kBB_;
    Teuchos::Array<UIntTuple> nodes_;
    Teuchos::Array<int> boundaryIndices_;
    Teuchos::Array<nodeType> nodeTypes_;
    Point origin_;

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
    nodeType
    getNodeType ( const direction prevDir,
                  const direction newDir
                ) const;

    direction
    getDirection ( const UIntTuple & node0,
                   const UIntTuple & node1
                 ) const;

    void
    updateGridDomainArea();

    Teuchos::RCP<Point>
    getX ( const UIntTuple & i ) const;

    bool
    boundaryStepper ( Teuchos::Array<UIntTuple>                  & boundaryNodes,
                      Teuchos::Array<direction>                  & directions,
                      const Teuchos::RCP<const Domain::Abstract> & domain
                    ) const;

    bool
    equal ( const UIntTuple & a,
            const UIntTuple & b
          ) const;

    UIntTuple
    findFirstBoundaryNode ( const Teuchos::RCP<const Domain::Abstract> & domain ) const;

    direction
    findNextDirection ( const IntTuple  & currentBoundaryNode,
                        const direction   prevDir
                      ) const;

    Teuchos::Tuple<direction,3>
    getNextDirections ( const direction dir ) const;

    UIntTuple
    step ( const UIntTuple                         & node,
           const direction                           dir,
           const Teuchos::RCP<const Domain::Abstract> & domain
         ) const;

    unsigned int
    i2kBoundingBox ( const UIntTuple & i ) const;

    Teuchos::RCP<UIntTuple>
    k2iBoundingBox ( const unsigned int k ) const;

    unsigned int
    kBoundingBox2kDomain ( const unsigned int k ) const;

    void
    pruneInitialTentacle ( Teuchos::Array<UIntTuple> & nodes,
                           Teuchos::Array<direction> & directions
                         ) const;
};

  } // namespace Grid
} // namespace Recti

#endif /* GRID_GENERAL_H_ */
