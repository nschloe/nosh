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
                                        unsigned int numGridPoints,
                                        unsigned int numBoundaryPoints ) :
GridVirtual( scaling,
             Teuchos::tuple<double>(h,h),
             gridDomainArea,
             numGridPoints,
             numBoundaryPoints )
{
}
// ============================================================================
GridUniformVirtual::GridUniformVirtual() :
GridVirtual()
{
}
// ============================================================================
GridUniformVirtual::~GridUniformVirtual()
{
}
// ============================================================================
double
GridUniformVirtual::getUniformH() const
{
  return h_[0];
}
// ============================================================================
