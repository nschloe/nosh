/*
 * GridReader.cpp
 *
 *  Created on: Nov 27, 2009
 *      Author: Nico Schl\"omer
 */

#include "GridUniformSquare.h"

#include "ioVirtual.h"
#include "ioFactory.h"

namespace GridReader {
// =============================================================================
void
read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
      const std::string                             & filePath,
            Teuchos::RCP<DoubleMultiVector>         & x,
            Teuchos::RCP<GridUniformVirtual>        & grid,
            Teuchos::ParameterList                  & params
    )
{
  // TODO Get some clues about which grid we read.
  //      Right now we can only read GridUniformSquare.

  Teuchos::RCP<GridUniformSquare> tmpGridUniformSquare = Teuchos::rcp( new GridUniformSquare() );
  tmpGridUniformSquare->read( Comm, filePath, x, params );
  grid = tmpGridUniformSquare; // slice
}
// =============================================================================
void
read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
      const std::string                             & filePath,
            Teuchos::RCP<ComplexMultiVector>        & x,
            Teuchos::RCP<GridUniformVirtual>        & grid,
            Teuchos::ParameterList                  & params
    )
{

  // TODO Get some clues about which grid we read.
  //      Right now we can only read GridUniformSquare.

  Teuchos::RCP<GridUniformSquare> tmpGridUniformSquare = Teuchos::rcp( new GridUniformSquare() );
  tmpGridUniformSquare->read( Comm, filePath, x, params );
  grid = tmpGridUniformSquare; // slice
}
// =============================================================================
} // namespace GridReader
