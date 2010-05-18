/*
 * GridReader.cpp
 *
 *  Created on: Nov 27, 2009
 *      Author: Nico Schl\"omer
 */

#include "Recti_Grid_Reader.h"

#include "Ginla_State.h"

#include "VIO_Reader_Factory.h"

// =============================================================================
void
Recti::Grid::Reader::
read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
      const std::string                             & filePath,
            Teuchos::RCP<Ginla::State>              & state,
            Teuchos::RCP<Recti::Grid::Uniform>      & grid,
            Teuchos::ParameterList                  & fieldData
    )
{
  Teuchos::RCP<VIO::Reader::Abstract> reader = VIO::Reader::Factory::create( filePath );
  
  // read all the values
  Teuchos::RCP<ComplexMultiVector> z;
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
  grid = Teuchos::rcp( new Uniform ( hh, numCells, bbIndex, boundaryIndices,
                                     scaling, origin ) );
                                     
  state = Teuchos::rcp( new Ginla::State( z, grid ) );

  return;
}
// =============================================================================
