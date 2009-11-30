/*
 * GridReader.cpp
 *
 *  Created on: Nov 27, 2009
 *      Author: Nico Schlšmer
 */

#include "GridReader.h"

#include "ioVirtual.h"
#include "ioFactory.h"

// =============================================================================
GridReader::GridReader()
{
}
// =============================================================================
GridReader::~GridReader()
{
}
// =============================================================================
void
GridReader::read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
                  const std::string                             & filePath,
                  Teuchos::RCP<DoubleMultiVector>               & x,
                  Teuchos::RCP<GridUniformVirtual>              & grid,
                  Teuchos::ParameterList                        & params
                ) const
{

  // TODO Get some clues about which grid we read.
  //      Right now we can only read GridSquare.

  Teuchos::RCP<GridSquare> tmpGridSquare = Teuchos::rcp( new GridSquare() );
  tmpGridSquare->read( Comm, filePath, x, params );
  grid = tmpGridSquare; // slice

}
// =============================================================================
