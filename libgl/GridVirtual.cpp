/*
 * GridVirtual.cpp
 *
 *  Created on: Nov 25, 2009
 *      Author: Nico Schlšmer
 */

#include "GridVirtual.h"

// ============================================================================
GridVirtual::GridVirtual( double scaling,
                          double h,
                          double gridDomainArea,
                          int    numGridPoints,
                          int    numBoundaryPoints ) :
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
  h_              *= ratio; // rescale h
  gridDomainArea_ *= ratio*ratio;  // rescale domain area
  scaling_         = scaling;
  return;
}
// ============================================================================
int
GridVirtual::getNumGridPoints() const
{
  return numGridPoints_;
}
// ============================================================================
int
GridVirtual::getNumBoundaryPoints() const
{
  return numBoundaryPoints_;
}
// ============================================================================
double
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
