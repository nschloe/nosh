/*
 * Grid.h
 *
 *  Created on: Nov 5, 2009
 *      Author: Nico Schl�mer
 */

#ifndef GRIDSQUARE_H_
#define GRIDSQUARE_H_

#include "GridVirtual.h"

#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Comm.hpp>

class GridSquare: virtual public GridVirtual
{

public:

    //! Default constructor.
    GridSquare ( const UIntTuple & numCells,
                 const double      edgeLength
               );

    //! Class constructor that only initializes the data members of this class.
    GridSquare ( const UIntTuple & numCells );

    virtual
    ~GridSquare();

    virtual Teuchos::RCP<DoubleTuple>
    getXLeft ( unsigned int k ) const; //!< Returns the value of \f$x\f$ left of point i.

    virtual Teuchos::RCP<DoubleTuple>
    getXRight ( unsigned int k ) const; //!< Returns the value of \f$x\f$ right of point i.

    virtual Teuchos::RCP<DoubleTuple>
    getXBelow ( unsigned int k ) const; //!< Returns the value of \f$x\f$ below point i.

    virtual Teuchos::RCP<DoubleTuple>
    getXAbove ( unsigned int k ) const; //!< Returns the value of \f$x\f$ above point i.


    virtual unsigned int
    getKLeft ( unsigned int k ) const; //!< Returns the running index \c k of the node left of \c i.

    virtual unsigned int
    getKRight ( unsigned int k ) const; //!< Returns the running index \c k of the node right of \c i.

    virtual unsigned int
    getKBelow ( unsigned int k ) const; //!< Returns the running index \c k of the node below \c i.

    virtual unsigned int
    getKAbove ( unsigned int k ) const; //!< Returns the running index \c k of the node above \c i.

    //! Returns the global index \ck of the \cl-th boundary node. Subsequent nodes \cl, \cl+1
    //! sit next to each other.
    virtual unsigned int
    boundaryIndex2globalIndex ( unsigned int l ) const;

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

    virtual void
    read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
           const std::string                             & filePath,
           Teuchos::RCP<DoubleMultiVector>               & x,
           Teuchos::ParameterList                        & params
         );

    virtual nodeType
    getNodeType ( unsigned int k ) const;

    virtual nodeType
    getBoundaryNodeType ( unsigned int l ) const;

    virtual double
    cellArea ( unsigned int k ) const;

protected:
    Teuchos::Tuple<unsigned int,2> numCells_; //!< Number of grid pieces in both x- and y-direction

    Teuchos::RCP<DoubleMultiVector>
    permuteLexicographic2Grid ( const DoubleMultiVector & xLexicographic
                              ) const;
    Teuchos::RCP<ComplexMultiVector>
    permuteLexicographic2Grid ( const ComplexMultiVector & xLexicographic
                              ) const;

    Teuchos::RCP<Epetra_MultiVector>
    permuteGrid2Lexicographic ( const Epetra_MultiVector & x
                              ) const;
    Teuchos::RCP<DoubleMultiVector>
    permuteGrid2Lexicographic ( const DoubleMultiVector & x
                              ) const;
    Teuchos::RCP<ComplexMultiVector>
    permuteGrid2Lexicographic ( const ComplexMultiVector & x
                              ) const;

private:

    //! Defines a subsequent order of boundary nodes by associating a running index \c l with
    //! a \f$d\f$-dimensional grid coordinate \f$i\f$.
    Teuchos::RCP<IntTuple>
    boundaryPosition ( unsigned int l ) const;

    Teuchos::RCP<DoubleTuple>
    getX ( IntTuple & i ) const;

    Teuchos::RCP<IntTuple>
    k2i ( unsigned int k ) const;

    //!< Converts a grid index i to a running index k
    unsigned int
    i2k ( Teuchos::RCP<IntTuple> & i ) const;

    nodeType
    getNodeTypeFromI ( IntTuple & i ) const;

private:
};

#endif /* GRIDSQUARE_H_ */
