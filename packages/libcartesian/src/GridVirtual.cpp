/*
 * GridVirtual.cpp
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schl\"omer
 */

#include "GridVirtual.h"

// ============================================================================
GridVirtual::GridVirtual ( const DoubleTuple  & h,
                           const double         gridDomainArea,
                           const unsigned int   numGridPoints,
                           const double         scaling ) :
        h_ ( h ),
        scaling_ ( scaling ),
        gridDomainArea_ ( gridDomainArea ),
        numGridPoints_ ( numGridPoints )
{
}
// ============================================================================
GridVirtual::GridVirtual() :
        h_ ( Teuchos::tuple<double> ( 0.0,0.0 ) ),
        scaling_ ( 1.0 ),
        gridDomainArea_ ( 0.0 ),
        numGridPoints_ ( 0 )
{
}
// ============================================================================
GridVirtual::~GridVirtual()
{
}
// ============================================================================
double
GridVirtual::getScaling() const
{
    return scaling_;
}
// ============================================================================
void
GridVirtual::updateScaling ( const double scaling )
{    
    TEUCHOS_ASSERT_INEQUALITY ( scaling, >, 0.0 );

    double ratio = scaling/scaling_;

    for ( int k=0; k<h_.size(); k++ )
        h_[k] *= ratio;
    gridDomainArea_ *= ratio*ratio;  // rescale domain area
    scaling_         = scaling;
    
    return;
}
// ============================================================================
DoubleTuple
GridVirtual::getH() const
{
    return h_;
}
// ============================================================================
unsigned int
GridVirtual::getNumGridPoints() const
{
    return numGridPoints_;
}
// ============================================================================
double
GridVirtual::getGridDomainArea() const
{
    return gridDomainArea_;
}
// =============================================================================
