/*
 * GridUniformSquare.h
 *
 *  Created on: Nov 30, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef GRID_UNIFORMSQUARE_H_
#define GRID_UNIFORMSQUARE_H_

#include "Recti_Grid_UniformAbstract.h"
#include "Recti_Grid_Square.h"

namespace Recti
{
  namespace Grid
  {

class UniformSquare:
            public UniformAbstract,
            public Square
{
public:
    //! Default constructor.
    UniformSquare ( const unsigned int numCells,
                    const double       edgeLength );

    virtual
    ~UniformSquare();

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

    void
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

    virtual void
    read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
           const std::string                             & filePath,
           Teuchos::RCP<ComplexMultiVector>              & x,
           Teuchos::ParameterList                        & params
         );

protected:
private:
};

  } // namespace Grid
} // namespace Recti

#endif /* GRID_UNIFORMSQUARE_H_ */
