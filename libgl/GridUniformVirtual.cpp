/*
 * GridUniformVirtual.cpp
 *
 *  Created on: Nov 30, 2009
 *      Author: Nico Schlšmer
 */

#include "GridUniformVirtual.h"

// ============================================================================
GridUniformVirtual::GridUniformVirtual( double scaling,
                                        double h,
                                        double gridDomainArea,
                                        int    numGridPoints,
                                        int    numBoundaryPoints ) :
GridVirtual( scaling,
             Teuchos::tuple(h,h),
             gridDomainArea,
             numGridPoints,
             numBoundaryPoints ),
h_( h )
{
}
// ============================================================================
GridUniformVirtual::~GridUniformVirtual()
{
}
// ============================================================================
void
GridUniformVirtual::setScaling( const double scaling )
{
  GridVirtual::setScaling( scaling );

  double ratio = scaling/scaling_;
  GridUniformVirtual::h_ *= ratio; // rescale h

  return;
}
// ============================================================================
double
GridUniformVirtual::getH() const
{
  return GridUniformVirtual::h_;
}
// ============================================================================
