/*
 * GridReader.cpp
 *
 *  Created on: Nov 27, 2009
 *      Author: Nico Schl\"omer
 */

#include "GridReader.h"

#include "ReaderFactory.h"

namespace GridReader {
// =============================================================================
void
read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
      const std::string                             & filePath,
            Teuchos::RCP<DoubleMultiVector>         & x,
            Teuchos::RCP<GridUniform>               & grid,
            Teuchos::ParameterList                  & params
    )
{
  TEST_FOR_EXCEPTION ( true,
                       std::logic_error,
                       "Not yet implemented." );

//   Teuchos::RCP<GridUniformSquare> tmpGridUniformSquare = Teuchos::rcp( new GridUniformSquare() );
//   tmpGridUniformSquare->read( Comm, filePath, x, params );
//   grid = tmpGridUniformSquare; // slice

//   Teuchos::RCP<VtiReader> reader = Teuchos::rcp( new VtiReader( filePath ) );
//   reader->read();
//   
//   return;
}
// =============================================================================
void
read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
      const std::string                             & filePath,
            Teuchos::RCP<ComplexMultiVector>        & z,
            Teuchos::RCP<GridUniform>               & grid,
            Teuchos::ParameterList                  & fieldData
    )
{
  
  Teuchos::RCP<AbstractImageReader> reader = ReaderFactory::create( filePath );
  
  // read all the values
  UIntTuple dims;
  DoubleTuple origin;
  DoubleTuple spacing; // see note below
  Teuchos::Array<int> bbIndex;
  reader->read( z, bbIndex, dims, origin, spacing, fieldData, Comm );
  
  // extract the necessary values
  double scaling = fieldData.get<double>( "scaling" );

  Teuchos::Array<int> boundaryIndices = fieldData.get<Teuchos::Array<int> >( "boundary indices" );
  
  // h and spacing should essentially deliver the same values, where SPACING is less reliable
  // as it may be stored with little precision in the file, depending on the file type.
  Teuchos::Array<double> hArray = fieldData.get<Teuchos::Array<double> >( "h" );

  TEUCHOS_ASSERT_EQUALITY( hArray.length(), 2 );
  DoubleTuple h = Teuchos::tuple( hArray[0], hArray[1] );
  
  TEST_FOR_EXCEPTION( fabs( h[0]-h[1] ) > 1.e-15,
                      std::logic_error,
                      "Spacing not uniform accross spatial dimensions: h = " << h << "." );

  // now create a grid out of what we got
  double hh = h[0];
  UIntTuple numCells = Teuchos::tuple(  dims[0]-1, dims[1]-1 );
  grid = Teuchos::rcp( new GridUniform ( hh, numCells, bbIndex, boundaryIndices,
                                         scaling, origin ) );

  return;
}
// =============================================================================
} // namespace GridReader
