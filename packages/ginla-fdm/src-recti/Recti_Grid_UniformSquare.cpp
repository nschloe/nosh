/*
 * GridUniformSquare.cpp
 *
 *  Created on: Nov 30, 2009
 *      Author: Nico Schl\"omer
 */

#include "Recti_Grid_UniformSquare.h"

// =============================================================================
// Class constructor
Recti::Grid::UniformSquare::
UniformSquare ( const unsigned int  numCells,
                const double        edgeLength ) :
        Abstract ( Teuchos::tuple ( edgeLength/numCells,
                                    edgeLength/numCells,
                                    0.0 ),
                   edgeLength * edgeLength,
                   (numCells+1) * (numCells+1)
                 ),
        UniformAbstract(),
        Square ( Teuchos::tuple( numCells,numCells ), edgeLength )
{
}
// =============================================================================
// Destructor
Recti::Grid::UniformSquare::
~UniformSquare()
{
}
// =============================================================================
void
Recti::Grid::UniformSquare::
writeWithGrid ( const Epetra_MultiVector     & x,
                const Teuchos::ParameterList & params,
                const std::string            & filePath
              ) const
{
  TEST_FOR_EXCEPTION( true,
                      std::logic_error,
                      "Not yet implemented." );

//     Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp ( IoFactory::createFileIo ( filePath ) );
// 
//     // append grid parameters
//     Teuchos::ParameterList extendedParams ( params );
//     extendedParams.get ( "scaling", scaling_ );
//     extendedParams.get ( "Nx", Nx_[0] );
// 
//     // reorder the grid to lexicographic ordering
//     Teuchos::RCP<Epetra_MultiVector> xLexicographic = permuteGrid2Lexicographic ( x );
// 
//     fileIo->write ( *xLexicographic, Nx_, GridSquare::h_, extendedParams );
}
// =============================================================================
void
Recti::Grid::UniformSquare::
writeWithGrid ( const DoubleMultiVector      & x,
                const Teuchos::ParameterList & params,
                const std::string            & filePath
              ) const
{
  TEST_FOR_EXCEPTION( true,
                      std::logic_error,
                      "Not yet implemented." );
                      
//     Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp ( IoFactory::createFileIo ( filePath ) );
// 
//     // append grid parameters
//     Teuchos::ParameterList extendedParams ( params );
//     extendedParams.get ( "scaling", scaling_ );
//     extendedParams.get ( "Nx", Nx_[0] );
//     // reorder the grid to lexicographic ordering
//     Teuchos::RCP<DoubleMultiVector> xLexicographic = permuteGrid2Lexicographic ( x );
// 
//     fileIo->write ( *xLexicographic, Nx_, GridSquare::h_, extendedParams );
}
// =============================================================================
void
Recti::Grid::UniformSquare::
writeWithGrid ( const ComplexMultiVector     & x,
                                   const Teuchos::ParameterList & params,
                                   const std::string            & filePath
                                 ) const
{
  TEST_FOR_EXCEPTION( true,
                      std::logic_error,
                      "Not yet implemented." );
  
//     Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp ( IoFactory::createFileIo ( filePath ) );
// 
//     // append grid parameters
//     Teuchos::ParameterList extendedParams ( params );
//     extendedParams.get ( "scaling", scaling_ );
//     extendedParams.get ( "Nx", Nx_[0] );
//     // reorder the grid to lexicographic ordering
//     Teuchos::RCP<ComplexMultiVector> xLexicographic = permuteGrid2Lexicographic ( x );
// 
//     fileIo->write ( *xLexicographic, Nx_, GridSquare::h_, extendedParams );
}
// =============================================================================
void
Recti::Grid::UniformSquare::
read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
                          const std::string                             & filePath,
                          Teuchos::RCP<DoubleMultiVector>               & x,
                          Teuchos::ParameterList                        & params
                        )
{
  TEST_FOR_EXCEPTION( true,
                      std::logic_error,
                      "Not yet implemented." );
                      
//     Teuchos::RCP<IoVirtual> fileIo = Teuchos::RCP<IoVirtual> (
//                                          IoFactory::createFileIo ( filePath ) );
// 
//     Teuchos::RCP<DoubleMultiVector> xLexicographic;
//     fileIo->read ( Comm, xLexicographic, params );
// 
//     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     // create the grid with the just attained information
//     TEST_FOR_EXCEPTION ( !params.isParameter ( "Nx" ),
//                          std::logic_error,
//                          "Parameter \"Nx\" not found in file " << filePath );
//     unsigned int numCells = params.get<int> ( "Nx" );
//     Nx_ = Teuchos::tuple<unsigned int> ( numCells, numCells );
// 
//     TEST_FOR_EXCEPTION ( !params.isParameter ( "scaling" ),
//                          std::logic_error,
//                          "Parameter \"scaling\" not found in file " << filePath );
// 
//     scaling_ = params.get<double> ( "scaling" );
// 
//     // initialization of the dependent members
//     double h = scaling_/numCells;
//     h_                 = Teuchos::tuple<double> ( h, h );
//     gridDomainArea_    = pow ( scaling_, 2 );
//     numGridPoints_     = ( Nx_[0]+1 ) * ( Nx_[1]+1 );
//     numBoundaryPoints_ = 2* ( Nx_[0]+Nx_[1] );
// 
//     // apply the grid ordering
//     x = permuteLexicographic2Grid ( *xLexicographic );
}
// =============================================================================
void
Recti::Grid::UniformSquare::
read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
                          const std::string                             & filePath,
                          Teuchos::RCP<ComplexMultiVector>                   & x,
                          Teuchos::ParameterList                        & params
                        )
{
  TEST_FOR_EXCEPTION( true,
                      std::logic_error,
                      "Not yet implemented." );

//     Teuchos::RCP<IoVirtual> fileIo = Teuchos::RCP<IoVirtual> (
//                                          IoFactory::createFileIo ( filePath ) );
// 
//     Teuchos::RCP<ComplexMultiVector> xLexicographic;
//     fileIo->read ( Comm, xLexicographic, params );
// 
//     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     // create the grid with the just attained information
//     TEST_FOR_EXCEPTION ( !params.isParameter ( "Nx" ),
//                          std::logic_error,
//                          "Parameter \"Nx\" not found in file " << filePath );
//     unsigned int numCells = params.get<int> ( "Nx" );
//     Nx_ = Teuchos::tuple<unsigned int> ( numCells, numCells );
// 
//     TEST_FOR_EXCEPTION ( !params.isParameter ( "scaling" ),
//                          std::logic_error,
//                          "Parameter \"scaling\" not found in file " << filePath );
// 
//     scaling_ = params.get<double> ( "scaling" );
// 
//     // initialization of the dependent members
//     double h = scaling_/numCells;
//     h_                 = Teuchos::tuple<double> ( h, h );
//     gridDomainArea_    = pow ( scaling_, 2 );
//     numGridPoints_     = ( Nx_[0]+1 ) * ( Nx_[1]+1 );
//     numBoundaryPoints_ = 2* ( Nx_[0]+Nx_[1] );
// 
//     // apply the grid ordering
//     x = permuteLexicographic2Grid ( *xLexicographic );
}
// =============================================================================
