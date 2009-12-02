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
scaling_(scaling),
h_( h ),
gridDomainArea_( gridDomainArea ),
numGridPoints_( numGridPoints ),
numBoundaryPoints_( numBoundaryPoints )
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

  double ratio = scaling/scaling_;
  for ( int k=0; k<h_.size(); k++)
    h_[k] *= ratio;
//  h_              *= ratio; // rescale h
  gridDomainArea_ *= ratio*ratio;  // rescale domain area
  scaling_         = scaling;
  return;
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
Teuchos::Tuple<double,2>
GridVirtual::getH() const
{
  return h_;
}
// ============================================================================
double
GridVirtual::getGridDomainArea() const
{
  return gridDomainArea_;
}
// =============================================================================
