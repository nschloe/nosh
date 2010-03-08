/*
 * Grid.cpp
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schl\"omer
 */

#include "Grid.h"

#include "VIO_Writer_Factory.h"

#include <Teuchos_ArrayView.hpp>

// ============================================================================
Grid::Grid ( const Teuchos::RCP<const DomainVirtual> & domain,
             const DoubleTuple                       & h
           ) :
        GridVirtual ( h, 0.0, 0 ),
        numCells_ ( Teuchos::tuple ( ( unsigned int ) 0, ( unsigned int ) 0 ) ),
        kBB_ ( Teuchos::Array<int>() ),
        nodes_ ( Teuchos::Array<UIntTuple>() ),
        boundaryIndices_ ( Teuchos::Array<int>() ),
        nodeTypes_ ( Teuchos::Array<nodeType>() ),
        origin_ ( Teuchos::tuple ( 0.0, 0.0 ) )
{
    // get the boundary box
    Teuchos::Tuple<double,4> bb = domain->getBoundingBox();

    origin_ = Teuchos::tuple ( bb[0], bb[1] );
    
    // lay a grid over the box
    numCells_[0] = ceil ( ( bb[2]-bb[0] ) / h_[0] );
    numCells_[1] = ceil ( ( bb[3]-bb[1] ) / h_[1] );
    
    unsigned int maxNumNodes = ( numCells_[0]+1 ) * ( numCells_[1]+1 );

    // find a boundary node by scanning through the bounding box, left to right,
    // bottom to top
    UIntTuple firstBoundaryNode = findFirstBoundaryNode ( domain );

    Teuchos::Array<direction> directions;
    directions.push_back ( RIGHT ); // arbitrarily chosen for bootstrapping
    nodes_.push_back ( firstBoundaryNode );

    boundaryStepper ( nodes_, directions, domain );

    // firstBoundaryNode possibly sat in a tentacle, and boundaryStepper
    // can't deal with that. Hence, manually remove this tentacle.
    pruneInitialTentacle ( nodes_, directions );

    unsigned int numBoundaryPoints = nodes_.length();

    boundaryIndices_.resize ( numBoundaryPoints );
    for ( unsigned int k=0; k<numBoundaryPoints; k++ )
        boundaryIndices_[k] = k;

    // get the node types out of the directions
    for ( unsigned int k=0; k<numBoundaryPoints-1; k++ )
        nodeTypes_.push_back ( getNodeType ( directions[k], directions[k+1] ) );

    // take care of the last element
    nodeTypes_.push_back ( getNodeType ( directions[numBoundaryPoints-1], directions[0] ) );

    // create connection between running bounding box indexing
    // and running domain indexing;
    // first, only for the boundary nodes, rest follows later on
    kBB_.resize ( maxNumNodes,-1 );
    for ( unsigned int k=0; k<numBoundaryPoints; k++ )
        kBB_[ i2kBoundingBox ( nodes_[k] ) ] = k;

    // Flood: Now, loop over the whole field, left to right, top to bottom, and check
    // the remaining nodes.
    // Watch out for tentacles:
    // Only admit those nodes which actually sit inside the boundary points.
    unsigned int k = nodes_.length(); // number of boundary nodes
    nodes_.resize ( maxNumNodes );
    nodeTypes_.resize ( maxNumNodes );
    UIntTuple node;
    int l = 0;

    Teuchos::Array<bool> leftSweep ( numCells_[0]+1 );
    Teuchos::Array<bool> rightSweep ( numCells_[0]+1 );
    Teuchos::Array<bool> isBoundary ( numCells_[0]+1 );

    // TODO Replace the following with a proper filling algorithm
    // http://en.wikipedia.org/wiki/Flood_fill
    for ( unsigned int j=0; j<numCells_[1]+1; j++ )
    {
        // fill isBoundary for the left/right sweeps
        for ( unsigned int i=0; i<numCells_[0]+1; i++ )
            isBoundary[i] = kBB_[l++]>=0;

        // sweep the row from left to right to find candidates for interior nodes
        bool maybeInside = false; // the left endpoint is *never inside (at most boundary)
        bool isSwitched = false;
        for ( unsigned int i=0; i<numCells_[0]+1; i++ )
        {
            // check if we are crossing the border
            if ( isBoundary[i] )
            {
                leftSweep[i] = false;
                if ( !isSwitched )
                {
                    maybeInside = !maybeInside;
                    isSwitched = true;
                }
                continue;
            }
            isSwitched = false; // re-set for the next switch
            leftSweep[i] = maybeInside;
        }

        // sweep the row from right to left to find candidates for interior nodes
        maybeInside = false; // the left endpoint is *never inside (at most boundary)
        isSwitched = false;
        for ( unsigned int i=0; i<numCells_[0]+1; i++ )
        {
            // check if we are crossing the border
            if ( isBoundary[numCells_[0]-i] )
            {
                rightSweep[numCells_[0]-i] = false;
                if ( !isSwitched )
                {
                    maybeInside = !maybeInside;
                    isSwitched = true;
                }
                continue;
            }
            isSwitched = false; // re-set for the next switch
            rightSweep[numCells_[0]-i] = maybeInside;
        }

        for ( unsigned int i=0; i<numCells_[0]+1; i++ )
        {
            // TODO checking for the domain should not be necessary.
            // but when running line by line, you can't tell otherwise
            if ( leftSweep[i] && rightSweep[i] )
            {
                node = Teuchos::tuple ( i,j );
                if ( domain->isInDomain ( *getX ( node ) ) )
                {
                    nodes_[k]     = node;
                    nodeTypes_[k] = GridVirtual::INTERIOR;
                    k++;
                }
            }
        }
    }

    numGridPoints_ = k;

    // all BOUNDARY_* and INTERIOR nodes have been classified now; cut off the rest
    nodes_.resize ( numGridPoints_ );
    nodeTypes_.resize ( numGridPoints_ );

    // update kBB with the interior nodes
    for ( unsigned int k=numBoundaryPoints; k<numGridPoints_; k++ )
        kBB_[ i2kBoundingBox ( nodes_[k] ) ] = k;

    updateGridDomainArea();

    return;
}
// ============================================================================
Grid::Grid ( const DoubleTuple         & h,
             const UIntTuple           & numCells,
             const Teuchos::Array<int> & kBB,
             const Teuchos::Array<int> & boundaryNodes,
             const double                scaling,
             const DoubleTuple         & origin
           ) :
        GridVirtual ( h, 0.0, 0, scaling ),
        numCells_ ( numCells ),
        kBB_ ( kBB ),
        nodes_ ( Teuchos::Array<UIntTuple>() ),
        boundaryIndices_ ( boundaryNodes ),
        nodeTypes_ ( Teuchos::Array<nodeType>() ),
        origin_ ( origin )
{ 
    int numBoundaryPoints = boundaryNodes.length();
  
    int numBoundaryBoxPoints = ( numCells_[0]+1 ) * ( numCells_[1]+1 );
    TEUCHOS_ASSERT_EQUALITY ( kBB_.length(), numBoundaryBoxPoints );

    // gather info about the i-j-location of the nodes
    nodes_.resize ( numBoundaryBoxPoints );
    numGridPoints_ = 0;
    unsigned int k = 0;
    for ( unsigned int j=0; j<numCells_[1]+1; j++ )
    {
        for ( unsigned int i=0; i<numCells_[0]+1; i++ )
        {
            if ( kBB_[k]>=0 )
            {
                numGridPoints_++;
                nodes_[kBB_[k]] = Teuchos::tuple ( i, j );
            }
            k++;
        }
    }
    nodes_.resize ( numGridPoints_ );

    // regenerate information about the node type
    direction firstDir = getDirection ( nodes_[boundaryIndices_[numBoundaryPoints-1]],
                                        nodes_[boundaryIndices_[0]] );
    direction currentDir;
    direction nextDir = firstDir;
    nodeTypes_ = Teuchos::Array<nodeType> ( numGridPoints_, INTERIOR ); // default: interior
    for ( int k=0; k<boundaryIndices_.length()-1; k++ )
    {
        currentDir = nextDir;
        nextDir = getDirection ( nodes_[boundaryIndices_[k]],
                                 nodes_[boundaryIndices_[k+1]] );
        nodeTypes_[k] = getNodeType ( currentDir, nextDir );
    }
    // set the type of the last node
    currentDir = nextDir;
    nextDir = firstDir;
    nodeTypes_[boundaryIndices_[numBoundaryPoints-1]] = getNodeType ( currentDir, nextDir );

    updateGridDomainArea();
    
    return;
}
// ============================================================================
Grid::~Grid()
{
}
// ============================================================================
void
Grid::updateGridDomainArea()
{
    gridDomainArea_ = 0.0;
    double h2 = h_[0]*h_[1];
    for ( int k=0; k<nodes_.length(); k++ )
    {
        switch ( nodeTypes_[k] )
        {
        case BOUNDARY_BOTTOMLEFTCONVEX:
        case BOUNDARY_TOPLEFTCONVEX:
        case BOUNDARY_BOTTOMRIGHTCONVEX:
        case BOUNDARY_TOPRIGHTCONVEX:
            gridDomainArea_ += 0.25 * h2;
            break;
        case BOUNDARY_BOTTOM:
        case BOUNDARY_TOP:
        case BOUNDARY_LEFT:
        case BOUNDARY_RIGHT:
            gridDomainArea_ += 0.5 * h2;
            break;
        case BOUNDARY_BOTTOMLEFTCONCAVE:
        case BOUNDARY_TOPLEFTCONCAVE:
        case BOUNDARY_BOTTOMRIGHTCONCAVE:
        case BOUNDARY_TOPRIGHTCONCAVE:
            gridDomainArea_ += 0.25 * h2;
            break;
        case INTERIOR:
            gridDomainArea_ += 1.0 * h2;
            break;
        default:
            TEST_FOR_EXCEPTION ( true,
                                 std::logic_error,
                                 "Illegal node type \"" << nodeTypes_[k] << "\"." );
        }
    }
}
// ============================================================================
Grid::direction
Grid::getDirection ( const UIntTuple & node0,
                     const UIntTuple & node1
                   ) const
{
    if ( node0[0] == node1[0] )
    {
        if ( node0[1]+1 == node1[1] )
            return UP;
        else if ( node0[1] == node1[1]+1 )
            return DOWN;
    }
    else if ( node0[1] == node1[1] )
    {
        if ( node0[0]+1 == node1[0] )
            return RIGHT;
        else if ( node0[0] == node1[0]+1 )
            return LEFT;
    }

    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "The nodes " << node0 << " and " << node1
                         << " do not sit next to each other." );
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getX( const unsigned int k ) const
{
    return getX( nodes_[k] );
}
// =============================================================================
Teuchos::RCP<DoubleTuple>
Grid::getX ( const UIntTuple & i ) const
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
Grid::getXLeft ( const UIntTuple & i ) const
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
Grid::getXRight ( const UIntTuple & i ) const
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
Grid::getXBelow ( const UIntTuple & i ) const
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
Grid::getXAbove ( const UIntTuple & i ) const
{
    Teuchos::RCP<DoubleTuple> x ( getX ( i ) );
    ( *x ) [1] += 0.5 * h_[1];
    return x;
}
// =============================================================================
unsigned int
Grid::getKLeft ( unsigned int kDomain ) const
{
    UIntTuple i = nodes_[kDomain];
    TEUCHOS_ASSERT_INEQUALITY( i[0], >, 0 );
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
    UIntTuple i = nodes_[kDomain];
    TEUCHOS_ASSERT_INEQUALITY( i[0], <, numCells_[0] );
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
    UIntTuple i = nodes_[kDomain];
    TEUCHOS_ASSERT_INEQUALITY( i[1], >, 0 );
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
    UIntTuple i = nodes_[kDomain];
    TEUCHOS_ASSERT_INEQUALITY( i[1], <, numCells_[1] );
    i[1]++;

    // get the running index of the bounding box, and translate it
    // to a running index in the domain
    int k = kBB_[ i2kBoundingBox ( i ) ];
    TEUCHOS_ASSERT_INEQUALITY ( k, >=, 0 );
    return k;
}
// =============================================================================
Teuchos::RCP<UIntTuple>
Grid::k2iBoundingBox ( const unsigned int k ) const
{
    return Teuchos::rcp ( new UIntTuple ( Teuchos::tuple<unsigned int> ( k%numCells_[0],k/numCells_[0] ) ) );
}
// =============================================================================
unsigned int
Grid::i2kBoundingBox ( const UIntTuple & i ) const
{
    return i[0] + i[1]* ( numCells_[0]+1 );
}
// =============================================================================
GridVirtual::nodeType
Grid::getNodeType ( unsigned int kDomain ) const
{
    return nodeTypes_[ kDomain ];
}
// =============================================================================
// TODO move to helpers, along with getX, getK, and so forth
UIntTuple
Grid::findFirstBoundaryNode ( const Teuchos::RCP<const DomainVirtual> & domain ) const
{  
    for ( unsigned int j=0; j<numCells_[1]; j++ )
        for ( unsigned int i=0; i<numCells_[0]; i++ )
            if ( domain->isInDomain ( *getX ( Teuchos::tuple ( i,j ) ) ) )
                return Teuchos::tuple ( i,j );

    // if you get here no node was found
    TEST_FOR_EXCEPTION ( true,
                         std::runtime_error,
                         "No boundary node was found."
                       );
}
// ============================================================================
bool
Grid::equal ( const UIntTuple & a,
              const UIntTuple & b
            ) const
{
    return a[0]==b[0] && a[1]==b[1];
}
// ============================================================================
bool
Grid::boundaryStepper ( Teuchos::Array<UIntTuple>               & boundaryNodes,
                        Teuchos::Array<direction>               & directions,
                        const Teuchos::RCP<const DomainVirtual> & domain
                      ) const
{
    // Try to step into direction, and return if not possible.
    UIntTuple nextNode;

    UIntTuple  node = boundaryNodes.back();
    direction prevDir = directions.back();

    direction newDir;

    Teuchos::Tuple<direction,3> nextDirections = getNextDirections ( prevDir );
    for ( int k=0; k<3; k++ ) // try going into any of the the three directions, in order
    {
        newDir = nextDirections[k];
        try
        {
            nextNode = step ( node, newDir, domain );
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
            if ( boundaryStepper ( boundaryNodes, directions, domain ) )
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
UIntTuple
Grid::step ( const UIntTuple                         & node,
             const direction                           dir,
             const Teuchos::RCP<const DomainVirtual> & domain
           ) const
{
    UIntTuple newNode ( node );
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
    if ( newNode[0]<0 || newNode[0]>numCells_[0]
         || newNode[1]<0 || newNode[1]>numCells_[1]
         ||  !domain->isInDomain ( *getX ( newNode ) ) )
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
Grid::writeWithGrid ( const Epetra_MultiVector     & x,
                      const Teuchos::ParameterList & params,
                      const std::string            & filePath
                    ) const
{
    Teuchos::RCP<VIO::Writer::Abstract> writer =
        VIO::Writer::Factory::createImageWriter ( filePath );

    // append grid parameters
    Teuchos::ParameterList extendedParams ( params );
    extendedParams.set ( "scaling", scaling_ );

    // Store h_ separately as one can not rely on the SPACING to be stored with
    // appropriate precision. Legacy VTK files, for example, reserve a mere seven,
    // eight significant digits for SPACING.
    Teuchos::Array<double> hArray( h_ );
    extendedParams.set ( "h" , hArray );
    
    writer->addParameterList ( extendedParams );
    writer->addFieldData ( boundaryIndices_, "boundaryNodes" );
    writer->setImageData ( x, numCells_, h_, kBB_ );
    writer->write();
}
// =============================================================================
void
Grid::writeWithGrid ( const DoubleMultiVector      & x,
                      const Teuchos::ParameterList & params,
                      const std::string            & filePath
                    ) const
{
    Teuchos::RCP<VIO::Writer::Abstract> writer =
        VIO::Writer::Factory::createImageWriter ( filePath );

    // append grid parameters
    Teuchos::ParameterList extendedParams ( params );
    extendedParams.set ( "scaling", scaling_ );
    
    // Store h_ separately as one can not rely on the SPACING to be stored with
    // appropriate precision. Legacy VTK files, for example, reserve a mere seven,
    // eight significant digits for SPACING.
    Teuchos::Array<double> hArray( h_ );
    extendedParams.set ( "h" , hArray );

    writer->addParameterList ( extendedParams );
    writer->addFieldData ( boundaryIndices_, "boundaryNodes" );
    writer->setImageData ( x, numCells_, h_, kBB_ );
    writer->write();
}
// =============================================================================
void
Grid::writeWithGrid ( const ComplexMultiVector     & z,
                      const Teuchos::ParameterList & params,
                      const std::string            & filePath
                    ) const
{
    Teuchos::RCP<VIO::Writer::Abstract> writer =
        VIO::Writer::Factory::createImageWriter ( filePath );

    // append grid parameters
    Teuchos::ParameterList extendedParams ( params );
    extendedParams.set ( "scaling", scaling_ );

    // Store h_ separately as one can not rely on the SPACING to be stored with
    // appropriate precision. Legacy VTK files, for example, reserve a mere seven,
    // eight significant digits for SPACING.
    Teuchos::Array<double> hArray( h_ );
    extendedParams.set ( "h" , hArray );
    
    writer->addParameterList ( extendedParams );
    writer->addFieldData ( boundaryIndices_, "boundary indices" );

    if ( z.getNumVectors() == 1 )
    {
        Teuchos::Array<std::string> zName ( 1, "psi" );
        writer->setImageData ( z, numCells_, h_, kBB_, zName );
    }
    else
    {
        writer->setImageData ( z, numCells_, h_, kBB_ );
    }
    writer->write();
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
const Teuchos::Array<int> &
Grid::getBoundaryIndices() const
{
   return boundaryIndices_;
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
Grid::pruneInitialTentacle ( Teuchos::Array<UIntTuple> & nodes,
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

    if ( lengthInitialTentacle>0 )
    {
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
    }

    return;
}
// =============================================================================
