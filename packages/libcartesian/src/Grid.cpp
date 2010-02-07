/*
 * Grid.cpp
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schl\"omer
 */

#include "Grid.h"

#include "ioVirtual.h"
#include "ioFactory.h"

#include <Teuchos_SerialDenseMatrix.hpp>

// ============================================================================
Grid::Grid ( const Teuchos::RCP<const DomainVirtual> & domain,
             const DoubleTuple                         h,
             const double                              scaling
           ) :
        domain_ ( domain ),
        h_ ( h ),
        scaling_ ( scaling )
{
    // get the boundary box
    Teuchos::Tuple<double,4> bb = domain_->getBoundingBox();

    origin_ = Teuchos::tuple ( bb[0], bb[1] );

    // lay a grid over the box
    Nx_[0] = ceil ( ( bb[2]-bb[0] ) / h_[0] );
    Nx_[1] = ceil ( ( bb[3]-bb[1] ) / h_[1] );

    unsigned int maxNumNodes = ( Nx_[0]+1 ) * ( Nx_[1]+1 );

    // find a boundary node by scanning through the bounding box, left to right,
    // bottom to top
    IntTuple firstBoundaryNode = findFirstBoundaryNode ();


    Teuchos::Array<direction> directions;
    directions.push_back ( RIGHT ); // arbitrarily chosen for bootstrapping
    nodes_.push_back ( firstBoundaryNode );
    boundaryStepper ( nodes_, directions );

    // firstBoundaryNode possibly sat in a tentacle, and boundaryStepper
    // can't deal with that. Hence, manually remove this tentacle.
    pruneInitialTentacle ( nodes_, directions );

    int numBoundaryPoints = nodes_.length();

    // get the node types out of the directions
    for ( int k=0; k<numBoundaryPoints-1; k++ )
        nodeTypes_.push_back ( getNodeType ( directions[k], directions[k+1] ) );

    // take care of the last element
    int k = numBoundaryPoints-1;
    nodeTypes_.push_back ( getNodeType ( directions[k], directions[0] ) );


    // create connection between running bounding box indexing
    // and running domain indexing;
    // first, only for the boundary nodes, rest follows later on
    kBB_.resize ( maxNumNodes,-1 );
    for ( k=0; k<numBoundaryPoints; k++ )
        kBB_[ i2kBoundingBox ( nodes_[k] ) ] = k;


    // Flood: Now, loop over the whole field, left to right, top to bottom, and check
    // the remaining nodes.
    // Watch out for tentacles:
    // Only admit those nodes which actually sit inside the boundary points.
    k = nodes_.length(); // number of boundary nodes
    nodes_.resize ( maxNumNodes );
    nodeTypes_.resize ( maxNumNodes );
    IntTuple node;
    int l = 0;
    for ( int j=0; j<Nx_[1]+1; j++ )
    {
        bool isInside = false; // the left endpoint is *never inside (at most boundary)
        bool isSwitched = false;
        for ( int i=0; i<Nx_[0]+1; i++ )
        {
            // check if we are crossing the border
            if ( kBB_[l++]>=0 )
            {
                if ( !isSwitched )
                {
                    isInside = !isInside;
                    isSwitched = true;
                }
                continue;
            }
            isSwitched = false; // re-set for the next switch

            node = Teuchos::tuple ( i,j );
            // TODO checking for the domain should not be necessary.
            // but when running line by line, you can't tell otherwise

            if ( isInside && domain_->isInDomain ( *getX ( node ) ) )
            {
                nodes_[k]     = node;
                nodeTypes_[k] = GridVirtual::INTERIOR;
                k++;
            }
        }
    }

    numGridPoints_ = k;

    // all BOUNDARY_* and INTERIOR nodes have been classified now; cut off the rest
    nodes_.resize ( numGridPoints_ );
    nodeTypes_.resize ( numGridPoints_ );

    // update kBB with the interior nodes
    for ( int k=numBoundaryPoints; k<numGridPoints_; k++ )
        kBB_[ i2kBoundingBox ( nodes_[k] ) ] = k;
}
// ============================================================================
Grid::Grid() :
        h_ ( Teuchos::tuple<double> ( 0.0,0.0 ) ),
        scaling_ ( 0.0 ),
        gridDomainArea_ ( 0.0 ),
        numGridPoints_ ( 0 ),
        numBoundaryPoints_ ( 0 )
{
}
// ============================================================================
Grid::~Grid()
{
}
// ============================================================================
double
Grid::getScaling() const
{
    return scaling_;
}
// ============================================================================
void
Grid::setScaling ( const double scaling )
{
    TEST_FOR_EXCEPTION ( scaling==0.0,
                         std::logic_error,
                         "Previous scaling value scaling_=0.0." );

    double ratio = scaling/scaling_;

    for ( int k=0; k<h_.size(); k++ )
        h_[k] *= ratio;
    gridDomainArea_ *= ratio*ratio;  // rescale domain area
    scaling_         = scaling;

    return;
}
// ============================================================================
DoubleTuple
Grid::getH() const
{
    return h_;
}
// ============================================================================
unsigned int
Grid::getNumGridPoints() const
{
    return numGridPoints_;
}
// ============================================================================
unsigned int
Grid::getNumBoundaryPoints() const
{
    return numBoundaryPoints_;
}
// ============================================================================
double
Grid::getGridDomainArea() const
{
    return gridDomainArea_;
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getX ( const IntTuple & i ) const
{
    Teuchos::RCP<DoubleTuple> x = Teuchos::rcp ( new DoubleTuple() );
    ( *x ) [0] = i[0] * h_[0] + origin_[0];
    ( *x ) [1] = i[1] * h_[1] + origin_[1];
    return x;
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getXLeft ( const unsigned int k ) const
{
    return getXLeft ( nodes_[k] );
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getXLeft ( const IntTuple & i ) const
{
    Teuchos::RCP<DoubleTuple> x ( getX ( i ) );
    ( *x ) [0] -= 0.5 * h_[0];
    return x;
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getXRight ( const unsigned int k ) const
{
    return getXRight ( nodes_[k] );
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getXRight ( const IntTuple & i ) const
{
    Teuchos::RCP<DoubleTuple> x ( getX ( i ) );
    ( *x ) [0] += 0.5 * h_[0];
    return x;
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getXBelow ( const unsigned int k ) const
{
    return getXBelow ( nodes_[k] );
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getXBelow ( const IntTuple & i ) const
{
    Teuchos::RCP<DoubleTuple> x ( getX ( i ) );
    ( *x ) [1] -= 0.5 * h_[1];
    return x;
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getXAbove ( const unsigned int k ) const
{
    return getXAbove ( nodes_[k] );
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getXAbove ( const IntTuple & i ) const
{
    Teuchos::RCP<DoubleTuple> x ( getX ( i ) );
    ( *x ) [1] += 0.5 * h_[1];
    return x;
}
// =============================================================================
unsigned int
Grid::getKLeft ( unsigned int kDomain ) const
{
    // get the left i
    IntTuple i = nodes_[kDomain];
    i[0]--;

    // get the running index of the bounding box, and translate it
    // to a running index in the domain
    int k = kBB_[ i2kBoundingBox ( i ) ];
    TEUCHOS_ASSERT_INEQUALITY ( k, >=, 0 );
    return k;
}
// =============================================================================
unsigned int
Grid::getKRight ( unsigned int kDomain ) const
{
    // get the left i
    IntTuple i = nodes_[kDomain];
    i[0]++;

    // get the running index of the bounding box, and translate it
    // to a running index in the domain
    int k = kBB_[ i2kBoundingBox ( i ) ];
    TEUCHOS_ASSERT_INEQUALITY ( k, >=, 0 );
    return k;
}
// =============================================================================
unsigned int
Grid::getKBelow ( unsigned int kDomain ) const
{
    // get the left i
    IntTuple i = nodes_[kDomain];
    i[1]--;

    // get the running index of the bounding box, and translate it
    // to a running index in the domain
    int k = kBB_[ i2kBoundingBox ( i ) ];
    TEUCHOS_ASSERT_INEQUALITY ( k, >=, 0 );
    return k;
}
// =============================================================================
unsigned int
Grid::getKAbove ( unsigned int kDomain ) const
{
    // get the left i
    IntTuple i = nodes_[kDomain];
    i[1]++;

    // get the running index of the bounding box, and translate it
    // to a running index in the domain
    int k = kBB_[ i2kBoundingBox ( i ) ];
    TEUCHOS_ASSERT_INEQUALITY ( k, >=, 0 );
    return k;
}
// =============================================================================
// maps a running index k to a 2D index i
// Teuchos::RCP<IntTuple >
// Grid::k2i ( unsigned int k ) const
// {
//     Teuchos::RCP<IntTuple > i = Teuchos::rcp ( new IntTuple() );
//
//     if ( k < Nx_[0] )
//     { // south
//         ( *i ) [0] = k;
//         ( *i ) [1] = 0;
//     }
//     else if ( k < Nx_[0]+Nx_[1] )
//     { // east
//         ( *i ) [0] = Nx_[0];
//         ( *i ) [1] = k - Nx_[0];
//     }
//     else if ( k < 2*Nx_[0]+Nx_[1] )
//     { // north
//         ( *i ) [0] = 2*Nx_[0]+Nx_[1] -k;
//         ( *i ) [1] = Nx_[1];
//     }
//     else if ( k < 2* ( Nx_[0]+Nx_[1] ) )
//     { // west
//         ( *i ) [0] = 0;
//         ( *i ) [1] = 2* ( Nx_[0]+Nx_[1] ) - k;
//     }
//     else if ( k < ( Nx_[0]+1 ) * ( Nx_[1]+1 ) )
//     { // on the interior
//         ( *i ) [0] = ( k - numBoundaryPoints_ ) % ( Nx_[0] - 1 ) + 1;
//         ( *i ) [1] = ( k - numBoundaryPoints_ ) / ( Nx_[0] - 1 ) + 1;
//     }
//     else
//     {
//         TEST_FOR_EXCEPTION ( true,
//                              std::logic_error,
//                              "Illegal running index   k = " << k );
//     }
//
//     return i;
// }
// =============================================================================
// maps a 2D index i to a running index k
// unsigned int
// Grid::i2k ( Teuchos::RCP<IntTuple> & i ) const
// {
//     int k;
//
//     if ( ( *i ) [1] == 0 ) // south
//         k = ( *i ) [0];
//     else if ( ( *i ) [0] == ( int ) Nx_[0] ) // east
//         k = ( *i ) [1] + Nx_[0];
//     else if ( ( *i ) [1] == ( int ) Nx_[1] ) // north
//         k = 2*Nx_[0]+Nx_[1] - ( *i ) [0];
//     else if ( ( *i ) [0] == 0 ) // west
//         k = 2* ( Nx_[0]+Nx_[1] ) - ( *i ) [1];
//     else if ( ( *i ) [0] > 0 && ( *i ) [0] < ( int ) Nx_[0] && ( *i ) [1] > 0 && ( *i ) [1] < ( int ) Nx_[1] ) // interior
//         k = 2* ( Nx_[0]+Nx_[1] )
//             + ( Nx_[0] - 1 ) * ( ( *i ) [1] - 1 )
//             + ( *i ) [0] - 1;
//     else
//         TEST_FOR_EXCEPTION ( true, std::logic_error,
//                              "Illegal 2D index   i = " << *i );
//
//     return k;
// }
// =============================================================================
Teuchos::RCP<IntTuple>
Grid::k2iBoundingBox ( const unsigned int k ) const
{
    return Teuchos::rcp ( new IntTuple ( Teuchos::tuple<int> ( k%Nx_[0],k/Nx_[0] ) ) );
}
// =============================================================================
unsigned int
Grid::i2kBoundingBox ( const IntTuple & i ) const
{
    return i[0] + i[1]* ( Nx_[0]+1 );
}
// =============================================================================
GridVirtual::nodeType
Grid::getNodeType ( unsigned int kDomain ) const
{
    return nodeTypes_[ kDomain ];
}
// =============================================================================
// GridVirtual::nodeType
// Grid::getNodeTypeFromI ( IntTuple & i ) const
// {
//     if ( i[0]==0 )
//     {
//         if ( i[1]==0 )
//             return GridVirtual::BOUNDARY_BOTTOMLEFTCONVEX;
//         else if ( i[1]== ( int ) Nx_[1] )
//             return GridVirtual::BOUNDARY_TOPLEFTCONVEX;
//         else
//             return GridVirtual::BOUNDARY_LEFT;
//     }
//     else if ( i[0]== ( int ) Nx_[1] )
//     {
//         if ( i[1]==0 )
//             return GridVirtual::BOUNDARY_BOTTOMRIGHTCONVEX;
//         else if ( i[1]== ( int ) Nx_[1] )
//             return GridVirtual::BOUNDARY_TOPRIGHTCONVEX;
//         else
//             return GridVirtual::BOUNDARY_RIGHT;
//     }
//     else if ( i[1]==0 )
//         return GridVirtual::BOUNDARY_BOTTOM;
//     else if ( i[1]== ( int ) Nx_[1] )
//         return GridVirtual::BOUNDARY_TOP;
//     else
//         return GridVirtual::INTERIOR;
// }
// ============================================================================
// TODO move to helpers, along with getX, getK, and so forth
IntTuple
Grid::findFirstBoundaryNode () const
{



    for ( unsigned int j=0; j<Nx_[1]; j++ )
        for ( unsigned int i=0; i<Nx_[0]; i++ )
            if ( domain_->isInDomain ( *getX ( Teuchos::tuple ( ( int ) i, ( int ) j ) ) ) )
                return Teuchos::tuple ( ( int ) i, ( int ) j );

    // if you get here no node was found
    TEST_FOR_EXCEPTION ( true,
                         std::runtime_error,
                         "No boundary node was found."
                       );
}
// ============================================================================
bool
Grid::equal ( const IntTuple & a,
              const IntTuple & b
            ) const
{
    return a[0]==b[0] && a[1]==b[1];
}
// ============================================================================
bool
Grid::boundaryStepper ( Teuchos::Array<IntTuple>  & boundaryNodes,
                        Teuchos::Array<direction> & directions
                      ) const
{
    // Try to step into direction, and return if not possible.
    IntTuple nextNode;

    IntTuple  & node    = boundaryNodes.back();
    direction & prevDir = directions.back();

    direction newDir;

    Teuchos::Tuple<direction,3> nextDirections = getNextDirections ( prevDir );
    for ( int k=0; k<3; k++ ) // try going into any of the the three directions, in order
    {
        newDir = nextDirections[k];
        try
        {
            nextNode = step ( node, newDir );
        }
        catch ( ... )
        {   // try next node
            continue;
        }

        // step successful
        if ( equal ( nextNode, boundaryNodes.front() ) )
        {   // the loop is closed
            directions[0] = newDir;
            return true;
        }
        else
        {   // append node to the list and restart
            boundaryNodes.push_back ( nextNode );
            directions.push_back ( newDir );
            if ( boundaryStepper ( boundaryNodes, directions ) )
                return true; // unwind
        }
    }

    // If getting here, all three steps were unsuccessful.
    // Pop the last node.
    boundaryNodes.pop_back();
    directions.pop_back();
    return false;
}
// ============================================================================
Teuchos::Tuple<Grid::direction,3>
Grid::getNextDirections ( const direction dir ) const
{
    switch ( dir )
    {
    case LEFT:
        return Teuchos::tuple ( UP, LEFT, DOWN );
    case RIGHT:
        return Teuchos::tuple ( DOWN, RIGHT, UP );
    case UP:
        return Teuchos::tuple ( RIGHT, UP, LEFT );
    case DOWN:
        return Teuchos::tuple ( LEFT, DOWN, RIGHT );
    default:
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Illegal direction \"" << dir << "\"." );
    }
}
// ============================================================================
// Grid::direction
// Grid::findNextDirection ( const IntTuple  & currentBoundaryNode,
//                           const direction   prevDir
//                         ) const
// {
//     Teuchos::RCP<DoubleTuple> leftPos  = getXLeft ( currentBoundaryNode );
//     Teuchos::RCP<DoubleTuple> rightPos = getXRight ( currentBoundaryNode );
//     Teuchos::RCP<DoubleTuple> abovePos = getXAbove ( currentBoundaryNode );
//     Teuchos::RCP<DoubleTuple> belowPos = getXBelow ( currentBoundaryNode );
//
//     switch ( prevDir )
//     {
//     case LEFT:
//         if ( currentBoundaryNode[1]!=Nx_[1] && domain_->isInDomain ( *abovePos ) )
//         {
//             return UP;
//         }
//         else if ( currentBoundaryNode[0]!=0 && domain_->isInDomain ( *leftPos ) )
//         {
//             return LEFT;
//         }
//         else if ( currentBoundaryNode[1]!=0 && domain_->isInDomain ( *belowPos ) )
//         {
//             return DOWN;
//         }
//         break;
//     case RIGHT:
//         std::cout << domain_->isInDomain ( *belowPos ) << std::endl;
//         std::cout << domain_->isInDomain ( *rightPos ) << std::endl;
//         std::cout << domain_->isInDomain ( *abovePos ) << std::endl;
//         if ( currentBoundaryNode[1]!=0 && domain_->isInDomain ( *belowPos ) )
//         {
//             return DOWN;
//         }
//         else if ( currentBoundaryNode[0]!=Nx_[0] && domain_->isInDomain ( *rightPos ) )
//         {
//             return RIGHT;
//         }
//         else if ( currentBoundaryNode[1]!=Nx_[1] && domain_->isInDomain ( *abovePos ) )
//         {
//             return UP;
//         }
//         break;
//     case UP:
//         if ( currentBoundaryNode[0]!=Nx_[0] && domain_->isInDomain ( *rightPos ) )
//         {
//             return RIGHT;
//         }
//         else if ( currentBoundaryNode[1]!=Nx_[1] && domain_->isInDomain ( *abovePos ) )
//         {
//             return UP;
//         }
//         else if ( currentBoundaryNode[0]!=0 && domain_->isInDomain ( *leftPos ) )
//         {
//             return LEFT;
//         }
//         break;
//     case DOWN:
//         if ( currentBoundaryNode[0]!=0 && domain_->isInDomain ( *leftPos ) )
//         {
//             return LEFT;
//         }
//         else if ( currentBoundaryNode[1]!=0 && domain_->isInDomain ( *belowPos ) )
//         {
//             return DOWN;
//         }
//         if ( currentBoundaryNode[0]!=Nx_[0] && domain_->isInDomain ( *rightPos ) )
//         {
//             return RIGHT;
//         }
//         break;
//     }
//
//     TEST_FOR_EXCEPTION ( true,
//                          std::logic_error,
//                          "Illegal direction \"" << prevDir << "\"." );
// }
// ============================================================================
IntTuple
Grid::step ( const IntTuple  & node,
             const direction   dir
           ) const
{
    IntTuple newNode ( node );
    switch ( dir )
    {
    case LEFT:
        newNode[0]--;
        break;
    case RIGHT:
        newNode[0]++;
        break;
    case DOWN:
        newNode[1]--;
        break;
    case UP:
        newNode[1]++;
        break;
    default:
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Illegal direction \"" << dir << "\"." );
    }

    // check of the new node sits in the domain, and throw exception if not
    if ( newNode[0]<0 || newNode[0]>Nx_[0]
         || newNode[1]<0 || newNode[1]>Nx_[1]
         ||  !domain_->isInDomain ( *getX ( newNode ) ) )
        throw std::exception();

    return newNode;
}
// ============================================================================
Grid::nodeType
Grid::getNodeType ( const direction prevDir,
                    const direction newDir
                  ) const
{
    switch ( prevDir )
    {
    case LEFT:
        switch ( newDir )
        {
        case LEFT:
            return BOUNDARY_TOP;
        case UP:
            return BOUNDARY_TOPRIGHTCONCAVE;
        case DOWN:
            return BOUNDARY_TOPLEFTCONVEX;
        case RIGHT:
            break;
        }
    case RIGHT:
        switch ( newDir )
        {
        case RIGHT:
            return BOUNDARY_BOTTOM;
        case UP:
            return BOUNDARY_BOTTOMRIGHTCONVEX;
        case DOWN:
            return BOUNDARY_BOTTOMLEFTCONCAVE;
        case LEFT:
            break;
        }
    case UP:
        switch ( newDir )
        {
        case UP:
            return BOUNDARY_RIGHT;
        case LEFT:
            return BOUNDARY_TOPRIGHTCONVEX;
        case RIGHT:
            return BOUNDARY_BOTTOMRIGHTCONCAVE;
        case DOWN:
            break;
        }
    case DOWN:
        switch ( newDir )
        {
        case DOWN:
            return BOUNDARY_LEFT;
        case LEFT:
            return BOUNDARY_TOPLEFTCONCAVE;
        case RIGHT:
            return BOUNDARY_BOTTOMLEFTCONVEX;
        case UP:
            break;
        }
    }

    // if you arrive here, something went wrong
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Illegal directions combination: previousDir="
                         << prevDir << ", nextDir=" << newDir << "." );
}
// =============================================================================
double
Grid::cellArea ( unsigned int k ) const
{
    switch ( nodeTypes_[k] )
    {
    case GridVirtual::INTERIOR:
        return h_[0]*h_[1];
        break;
    case GridVirtual::BOUNDARY_BOTTOM:
    case GridVirtual::BOUNDARY_TOP:
    case GridVirtual::BOUNDARY_LEFT:
    case GridVirtual::BOUNDARY_RIGHT:
        return 0.5*h_[0]*h_[1];
        break;
    case GridVirtual::BOUNDARY_BOTTOMLEFTCONCAVE:
    case GridVirtual::BOUNDARY_BOTTOMRIGHTCONCAVE:
    case GridVirtual::BOUNDARY_TOPLEFTCONCAVE:
    case GridVirtual::BOUNDARY_TOPRIGHTCONCAVE:
        return 0.75*h_[0]*h_[1];
    case GridVirtual::BOUNDARY_BOTTOMLEFTCONVEX:
    case GridVirtual::BOUNDARY_BOTTOMRIGHTCONVEX:
    case GridVirtual::BOUNDARY_TOPLEFTCONVEX:
    case GridVirtual::BOUNDARY_TOPRIGHTCONVEX:
        return 0.25*h_[0]*h_[1];
        break;
    }

    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Illegal node type \""<<  nodeTypes_[k] << "\"." );
}
// =============================================================================
void
Grid::read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
             const std::string                             & filePath,
             Teuchos::RCP<DoubleMultiVector>               & x,
             Teuchos::ParameterList                        & params
           )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
void
Grid::read ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
             const std::string                             & filePath,
             Teuchos::RCP<ComplexMultiVector>              & x,
             Teuchos::ParameterList                        & params
           )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
void
Grid::writeWithGrid ( const Epetra_MultiVector     & x,
                      const Teuchos::ParameterList & params,
                      const std::string            & filePath
                    ) const
{
    Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp ( IoFactory::createFileIo ( filePath ) );

    // append grid parameters
    Teuchos::ParameterList extendedParams ( params );
    extendedParams.get ( "scaling", scaling_ );

    // create a list with the x,y values of the points
    Teuchos::Array<DoubleTuple> loc ( numGridPoints_ );
    for ( int k; k<numGridPoints_; k++ )
        loc[k] = *getX ( nodes_[k] );

    fileIo->write ( x, loc, extendedParams );
}
// =============================================================================
void
Grid::writeWithBoundingBoxGrid ( const Epetra_MultiVector     & x,
                                 const Teuchos::ParameterList & params,
                                 const std::string            & filePath
                               ) const
{
    Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp ( IoFactory::createFileIo ( filePath ) );

    // append grid parameters
    Teuchos::ParameterList extendedParams ( params );
    extendedParams.get ( "scaling", scaling_ );

    double dummyValue = 0.0;
    fileIo->write ( x, Nx_, h_, kBB_, extendedParams, dummyValue );
}
// =============================================================================
void
Grid::writeWithGrid ( const DoubleMultiVector      & x,
                      const Teuchos::ParameterList & params,
                      const std::string            & filePath
                    ) const
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
void
Grid::writeWithGrid ( const ComplexMultiVector     & x,
                      const Teuchos::ParameterList & params,
                      const std::string            & filePath
                    ) const
{
    Teuchos::RCP<IoVirtual> fileIo = Teuchos::rcp ( IoFactory::createFileIo ( filePath ) );

    // append grid parameters
    Teuchos::ParameterList extendedParams ( params );
    extendedParams.get ( "scaling", scaling_ );

    double dummyValue = 0.0;
    fileIo->write ( x, Nx_, h_, kBB_, extendedParams, dummyValue );
}
// =============================================================================
GridVirtual::nodeType
Grid::getBoundaryNodeType ( unsigned int l ) const
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
unsigned int
Grid::boundaryIndex2globalIndex ( unsigned int l ) const
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
void
Grid::pruneInitialTentacle ( Teuchos::Array<IntTuple>  & nodes,
                             Teuchos::Array<direction> & directions
                           ) const
{
    int numNodes = nodes.length();
    int lengthInitialTentacle = 0;
    for ( int k=1; k<numNodes; k++ )
    {
        if ( equal ( nodes[k], nodes[numNodes-k] ) )
            lengthInitialTentacle++;
        else
            break;
    }
    
    // Delete the respective nodes at the beginning and the end.
    // Note that one more node is deleted at the beginning (the first node),
    // *and* one more node is preserved at the end, the tentacle root.
    for ( int k=0; k<lengthInitialTentacle-1; k++ )
    {
        nodes.pop_back();
        directions.pop_back();
    }
    for ( int k=0; k<lengthInitialTentacle+1; k++ )
    {
        nodes.erase ( nodes.begin() );
        directions.erase ( directions.begin() );
    }

    return;
}
// =============================================================================
