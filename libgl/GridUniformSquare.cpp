/*
 * GridUniformSquare.cpp
 *
 *  Created on: Nov 30, 2009
 *      Author: Nico Schlšmer
 */

#include "GridUniformSquare.h"

#include "ioVirtual.h"
#include "ioFactory.h"

// =============================================================================
// Class constructor
GridUniformSquare::GridUniformSquare( unsigned int nx, double scaling) :
  GridVirtual( scaling,
               Teuchos::tuple<double>( scaling/nx,scaling/nx ),
               pow( scaling, 2 ),
               pow( nx+1, 2 ),
               4*nx
             ),
  GridUniformVirtual(),
  GridSquare( Teuchos::tuple<unsigned int>(nx,nx) )
{
}
// =============================================================================
// Destructor
GridUniformSquare::~GridUniformSquare()
{
}
// =============================================================================
void
GridUniformSquare::writeWithGrid( const DoubleMultiVector      & x,
                                  const Teuchos::ParameterList & params,
                                  const std::string            & filePath
                                ) const
{
  Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp(IoFactory::createFileIo(filePath));

  // append grid parameters
  Teuchos::ParameterList extendedParams( params );
  extendedParams.get("scaling", scaling_ );
  extendedParams.get("Nx", Nx_[0] );

  // reorder the grid to lexicographic ordering
  Teuchos::RCP<DoubleMultiVector> xLexicographic = permuteGrid2Lexicographic(x);

  fileIo->write( *xLexicographic, Nx_, GridSquare::h_, extendedParams);
}
// =============================================================================
void
GridUniformSquare::read( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
                         const std::string                             & filePath,
                         Teuchos::RCP<DoubleMultiVector>               & x,
                         Teuchos::ParameterList                        & params
                       )
{
  Teuchos::RCP<IoVirtual> fileIo = Teuchos::RCP<IoVirtual>(
                    IoFactory::createFileIo(filePath));

  Teuchos::RCP<DoubleMultiVector> xLexicographic;
  fileIo->read(Comm, xLexicographic, params);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // create the grid with the just attained information
  TEST_FOR_EXCEPTION( !params.isParameter("Nx"),
                      std::logic_error,
                      "Parameter \"Nx\" not found in file " << filePath );
  unsigned int nx = params.get<int>("Nx");
  Nx_ = Teuchos::tuple<unsigned int>( nx, nx );

  TEST_FOR_EXCEPTION( !params.isParameter("scaling"),
                      std::logic_error,
                      "Parameter \"scaling\" not found in file " << filePath );

  scaling_ = params.get<double>("scaling");

  // initialization of the dependent members
  double h = scaling_/nx;
  h_                 = Teuchos::tuple<double>( h, h );
  gridDomainArea_    = pow( scaling_, 2 );
  numGridPoints_     = (Nx_[0]+1)*(Nx_[1]+1);
  numBoundaryPoints_ = 2*(Nx_[0]+Nx_[1]);

  // apply the grid ordering
  x = permuteLexicographic2Grid( *xLexicographic );
}
// =============================================================================
