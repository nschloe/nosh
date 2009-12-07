/*
 * GridVirtual.cpp
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schlšmer
 */

#include "GridVirtual.h"

// ============================================================================
GridVirtual::GridVirtual( double scaling,
                          Teuchos::Tuple<double,2> h,
                          double gridDomainArea,
                          unsigned int numGridPoints,
                          unsigned int numBoundaryPoints ) :
h_( h ),
scaling_( scaling ),
gridDomainArea_( gridDomainArea ),
numGridPoints_( numGridPoints ),
numBoundaryPoints_( numBoundaryPoints )
{
}
// ============================================================================
GridVirtual::GridVirtual() :
h_( Teuchos::tuple<double>(0.0,0.0) ),
scaling_( 0.0 ),
gridDomainArea_( 0.0 ),
numGridPoints_( 0 ),
numBoundaryPoints_( 0 )
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
GridVirtual::setScaling( const double scaling )
{
  TEST_FOR_EXCEPTION( scaling==0.0,
                      std::logic_error,
                      "Previous scaling value scaling_=0.0." );

  double ratio = scaling/scaling_;

  for ( int k=0; k<h_.size(); k++)
    h_[k] *= ratio;
  gridDomainArea_ *= ratio*ratio;  // rescale domain area
  scaling_         = scaling;

  return;
}
// ============================================================================
Teuchos::Tuple<double,2>
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
unsigned int
GridVirtual::getNumBoundaryPoints() const
{
  return numBoundaryPoints_;
}
// ============================================================================
double
GridVirtual::getGridDomainArea() const
{
  return gridDomainArea_;
}
// =============================================================================
