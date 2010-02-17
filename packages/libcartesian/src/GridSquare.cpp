/*
 * Grid.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: Nico Schl\"omer
 */
#include <Epetra_Vector.h>

#include "GridSquare.h"

// =============================================================================
// Class constructor
GridSquare::GridSquare ( const UIntTuple & numCells,
                         const double      edgeLength ) :
        GridVirtual ( Teuchos::tuple<double> ( edgeLength/numCells[0],edgeLength/numCells[1] ),
                      pow ( edgeLength, 2 ),
                      ( numCells[0]+1 ) * ( numCells[1]+1 )
                    ),
        numCells_ ( numCells )
{
}
// =============================================================================
// Class constructor for this classes data only.
GridSquare::GridSquare ( const UIntTuple & numCells ) :
        GridVirtual(),
        numCells_ ( numCells )
{
}
// =============================================================================
// Destructor
GridSquare::~GridSquare()
{
}
// =============================================================================
Teuchos::RCP<IntTuple>
GridSquare::boundaryPosition ( unsigned int l ) const
{
    Teuchos::RCP<IntTuple> i = Teuchos::rcp ( new IntTuple() );

    unsigned int numBoundaryPoints = 2 * ( numCells_[0] + numCells_[1] );

    // start at the bottom left, and go around counter-clockwise
    if ( l < numCells_[0] )
    { // south
        ( *i ) [0] = l;
        ( *i ) [1] = 0;
    }
    else if ( l < numCells_[0]+numCells_[1] )
    { // east
        ( *i ) [0] = numCells_[0];
        ( *i ) [1] = l - numCells_[0];
    }
    else if ( l < 2*numCells_[0]+numCells_[1] )
    { // north
        ( *i ) [0] = 2*numCells_[0]+numCells_[1] - l;
        ( *i ) [1] = numCells_[1];
    }
    else if ( l < numBoundaryPoints )
    { // west
        ( *i ) [0] = 0;
        ( *i ) [1] = 2* ( numCells_[0]+numCells_[1] )  - l;
    }
    else
    {
        TEST_FOR_EXCEPTION ( true, std::logic_error,
                             "Given index l=" << l
                             << "larger than the number of boundary nodes n="
                             << numBoundaryPoints  << "." );
    }
    return i;
}
// =============================================================================
unsigned int
GridSquare::boundaryIndex2globalIndex ( unsigned int l ) const
{
    Teuchos::RCP<IntTuple> i = boundaryPosition ( l );
    return i2k ( i );
}
// =============================================================================
GridSquare::nodeType
GridSquare::getNodeType ( unsigned int k ) const
{
    return getNodeTypeFromI ( * ( k2i ( k ) ) );
}
// =============================================================================
GridSquare::nodeType
GridSquare::getBoundaryNodeType ( unsigned int l ) const
{
    return getNodeTypeFromI ( * ( boundaryPosition ( l ) ) );
}
// =============================================================================
double
GridSquare::cellArea ( unsigned int k ) const
{
    unsigned int numBoundaryPoints = 2 * ( numCells_[0] + numCells_[1] );
  
    // TODO Use nodeType association here and possibly move out to Virtual
    if ( k == 0 || k == numCells_[0] || k == numCells_[0]+numCells_[1] || k == 2*numCells_[0]+numCells_[1] )
        // corner
        return 0.25*h_[0]*h_[1];
    else if ( k < numBoundaryPoints )
        // edge
        return 0.5*h_[0]*h_[1];
    else
        // interior
        return h_[0]*h_[1];
}
// =============================================================================
// TODO implement this using import/export mechanisms
Teuchos::RCP<Epetra_MultiVector>
GridSquare::permuteGrid2Lexicographic ( const Epetra_MultiVector & x
                                      ) const
{
    TEUCHOS_ASSERT_EQUALITY ( ( unsigned int ) x.GlobalLength(), numGridPoints_ );

    int numVectors = x.NumVectors();
    Teuchos::RCP<Epetra_MultiVector> xLexicographic =
        Teuchos::rcp ( new Epetra_MultiVector ( x.Map(),numVectors ) );

    // Loop through the lexicographic ordering.
    for ( int l=0; l<numVectors; l++ )
    {
        int k = 0;
        Teuchos::RCP<IntTuple> index = Teuchos::rcp ( new IntTuple() );
        for ( unsigned int j = 0; j < numCells_[1] + 1; j++ )
        {
            ( *index ) [1] = j;
            for ( unsigned int i = 0; i < numCells_[0] + 1; i++ )
            {
                ( *index ) [0] = i;
                int kGlobal = i2k ( index );
                xLexicographic->ReplaceGlobalValue ( k++, l, ( *x ( l ) ) [kGlobal] );
            }
        }
    }

    return xLexicographic;
}
// =============================================================================
// TODO implement this using import/export mechanisms
// TODO templatetize this
Teuchos::RCP<DoubleMultiVector>
GridSquare::permuteGrid2Lexicographic ( const DoubleMultiVector & x
                                      ) const
{
    TEUCHOS_ASSERT_EQUALITY ( x.getGlobalLength(), numGridPoints_ );

    unsigned int numVectors = x.getNumVectors();
    Teuchos::RCP<DoubleMultiVector> xLexicographic =
        Teuchos::rcp ( new DoubleMultiVector ( x.getMap(),numVectors ) );

    // Loop through the lexicographic ordering.
    for ( unsigned int l=0; l<numVectors; l++ )
    {
        Teuchos::ArrayRCP<const double> xView = x.getVector ( l )->get1dView();
        int k = 0;
        Teuchos::RCP<IntTuple> index = Teuchos::rcp ( new IntTuple() );
        for ( unsigned int j = 0; j < numCells_[1] + 1; j++ )
        {
            ( *index ) [1] = j;
            for ( unsigned int i = 0; i < numCells_[0] + 1; i++ )
            {
                ( *index ) [0] = i;
                int kGlobal = i2k ( index );
                xLexicographic->replaceGlobalValue ( k++, l, xView[kGlobal] );
            }
        }
    }

    return xLexicographic;
}
// =============================================================================
// TODO implement this using import/export mechanisms
// TODO templatetize this
Teuchos::RCP<ComplexMultiVector>
GridSquare::permuteGrid2Lexicographic ( const ComplexMultiVector & x
                                      ) const
{
    TEUCHOS_ASSERT_EQUALITY ( x.getGlobalLength(), numGridPoints_ );

    unsigned int numVectors = x.getNumVectors();
    Teuchos::RCP<ComplexMultiVector> xLexicographic =
        Teuchos::rcp ( new ComplexMultiVector ( x.getMap(),numVectors ) );

    // Loop through the lexicographic ordering.
    for ( unsigned int l=0; l<numVectors; l++ )
    {
        Teuchos::ArrayRCP<const std::complex<double> > xView = x.getVector ( l )->get1dView();
        int k = 0;
        Teuchos::RCP<IntTuple> index = Teuchos::rcp ( new IntTuple() );
        for ( unsigned int j = 0; j < numCells_[1] + 1; j++ )
        {
            ( *index ) [1] = j;
            for ( unsigned int i = 0; i < numCells_[0] + 1; i++ )
            {
                ( *index ) [0] = i;
                int kGlobal = i2k ( index );
                xLexicographic->replaceGlobalValue ( k++, l, xView[kGlobal] );
            }
        }
    }

    return xLexicographic;
}
// =============================================================================
//// TODO implement this using import/export mechanisms
//// TODO templatetize this
Teuchos::RCP<DoubleMultiVector>
GridSquare::permuteLexicographic2Grid ( const DoubleMultiVector & xLexicographic
                                      ) const
{
    TEST_FOR_EXCEPTION ( xLexicographic.getGlobalLength() != numGridPoints_,
                         std::logic_error,
                         "Global length of the input vector x ("
                         << xLexicographic.getGlobalLength() << ") does not coincide "
                         << "with with number of unknowns on the grid ("
                         << numGridPoints_ << ")." );

    unsigned int numVectors = xLexicographic.getNumVectors();
    Teuchos::RCP<DoubleMultiVector> x =
        Teuchos::rcp ( new DoubleMultiVector ( xLexicographic.getMap(),numVectors ) );

    // Loop through the lexicographic ordering.
    for ( unsigned int l=0; l<numVectors; l++ )
    {
        Teuchos::ArrayRCP<const double> xLexicographicView = xLexicographic.getVector ( l )->get1dView();
        int k = 0;
        Teuchos::RCP<IntTuple> index = Teuchos::rcp ( new IntTuple() );
        for ( unsigned int j = 0; j < numCells_[1] + 1; j++ )
        {
            ( *index ) [1] = j;
            for ( unsigned int i = 0; i < numCells_[0] + 1; i++ )
            {
                ( *index ) [0] = i;
                int kGlobal = i2k ( index );
                x->replaceGlobalValue ( kGlobal, l, xLexicographicView[k++] );
            }
        }
    }

    return x;
}
// =============================================================================
//// TODO implement this using import/export mechanisms
//// TODO templatetize this
Teuchos::RCP<ComplexMultiVector>
GridSquare::permuteLexicographic2Grid ( const ComplexMultiVector & xLexicographic
                                      ) const
{
    TEST_FOR_EXCEPTION ( xLexicographic.getGlobalLength() != numGridPoints_,
                         std::logic_error,
                         "Global length of the input vector x ("
                         << xLexicographic.getGlobalLength() << ") does not coincide "
                         << "with with number of unknowns on the grid ("
                         << numGridPoints_ << ")." );

    unsigned int numVectors = xLexicographic.getNumVectors();
    Teuchos::RCP<ComplexMultiVector> x =
        Teuchos::rcp ( new ComplexMultiVector ( xLexicographic.getMap(),numVectors ) );

    // Loop through the lexicographic ordering.
    for ( unsigned int l=0; l<numVectors; l++ )
    {
        Teuchos::ArrayRCP<const std::complex<double> > xLexicographicView = xLexicographic.getVector ( l )->get1dView();
        int k = 0;
        Teuchos::RCP<IntTuple> index = Teuchos::rcp ( new IntTuple() );
        for ( unsigned int j = 0; j < numCells_[1] + 1; j++ )
        {
            ( *index ) [1] = j;
            for ( unsigned int i = 0; i < numCells_[0] + 1; i++ )
            {
                ( *index ) [0] = i;
                int kGlobal = i2k ( index );
                x->replaceGlobalValue ( kGlobal, l, xLexicographicView[k++] );
            }
        }
    }

    return x;
}
// =============================================================================
void
GridSquare::writeWithGrid ( const DoubleMultiVector      & x,
                            const Teuchos::ParameterList & params,
                            const std::string            & filePath
                          ) const
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );

//     Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp ( IoFactory::createFileIo ( filePath ) );
//
//     // append grid parameters
//     Teuchos::ParameterList extendedParams ( params );
//     extendedParams.get ( "scaling", scaling_ );
//     extendedParams.get ( "numCells", numCells_[0] );
//     extendedParams.get ( "Ny", numCells_[1] );
//
//     // reorder the grid to lexicographic ordering
//     Teuchos::RCP<DoubleMultiVector> xLexicographic = permuteGrid2Lexicographic ( x );
//
//     fileIo->write ( *xLexicographic, numCells_, h_, extendedParams );
}
// =============================================================================
void
GridSquare::writeWithGrid ( const ComplexMultiVector     & x,
                            const Teuchos::ParameterList & params,
                            const std::string            & filePath
                          ) const
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );

//     Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp ( IoFactory::createFileIo ( filePath ) );
//
//     // append grid parameters
//     Teuchos::ParameterList extendedParams ( params );
//     extendedParams.get ( "scaling", scaling_ );
//     extendedParams.get ( "numCells", numCells_[0] );
//     extendedParams.get ( "Ny", numCells_[1] );
//
//     // reorder the grid to lexicographic ordering
//     Teuchos::RCP<ComplexMultiVector> xLexicographic = permuteGrid2Lexicographic ( x );
//
//     fileIo->write ( *xLexicographic, numCells_, h_, extendedParams );
}
// =============================================================================
void
GridSquare::read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
                   const std::string                             & filePath,
                   Teuchos::RCP<DoubleMultiVector>               & x,
                   Teuchos::ParameterList                        & params
                 )
{
    TEST_FOR_EXCEPTION ( true,
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
//     TEST_FOR_EXCEPTION ( !params.isParameter ( "numCells" ),
//                          std::logic_error,
//                          "Parameter \"numCells\" not found." );
//     TEST_FOR_EXCEPTION ( !params.isParameter ( "Ny" ),
//                          std::logic_error,
//                          "Parameter \"Ny\" not found." );
//     numCells_ = Teuchos::tuple<unsigned int> ( params.get<unsigned int> ( "numCells" ), params.get<unsigned int> ( "Ny" ) );
//
//
//     TEST_FOR_EXCEPTION ( !params.isParameter ( "scaling" ),
//                          std::logic_error,
//                          "Parameter \"scaling\" not found." );
//     scaling_ = params.get<double> ( "scaling" );
//
//     // initialization of the dependent members
//     h_                 = Teuchos::tuple<double> ( scaling_/numCells_[0], scaling_/numCells_[1] );
//     gridDomainArea_    = pow ( scaling_, 2 );
//     numGridPoints_     = ( numCells_[0]+1 ) * ( numCells_[1]+1 );
//     numBoundaryPoints_ = 2* ( numCells_[0]+numCells_[1] );
//
//     TEST_FOR_EXCEPTION ( !x.is_valid_ptr() || x.is_null(),
//                          std::logic_error,
//                          "x not properly initialized." );
//
//     // apply the grid ordering
//     x = permuteLexicographic2Grid ( *xLexicographic );
}
// =============================================================================
