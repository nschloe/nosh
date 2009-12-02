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
  GridUniformVirtual( scaling,
                      scaling / nx,
                      pow( scaling, 2 ),
                     (nx+1)*(nx+1), 4*nx
                    ),
  GridSquare( Teuchos::tuple<unsigned int>(nx,nx),
              scaling )
{
}
// =============================================================================
// Destructor
GridUniformSquare::~GridUniformSquare()
{
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
                      "Parameter \"Nx\" not found." );
  unsigned int nx = params.get<int>("Nx");
  Nx_ = Teuchos::tuple<unsigned int>( nx, nx );

  TEST_FOR_EXCEPTION( !params.isParameter("scaling"),
                      std::logic_error,
                      "Parameter \"scaling\" not found." );

  scaling_ = params.get<double>("scaling");

  // initialization of the dependent members
  double h = scaling_/nx;
  GridUniformVirtual::h_ = h;
  GridSquare::h_         = Teuchos::tuple<double>( h, h );
  gridDomainArea_    = pow( scaling_, 2 );
  numGridPoints_     = (Nx_[0]+1)*(Nx_[1]+1);
  numBoundaryPoints_ = 2*(Nx_[0]+Nx_[1]);

  // apply the grid ordering
  x = permuteLexicographic2Grid( *xLexicographic );
}
// =============================================================================
