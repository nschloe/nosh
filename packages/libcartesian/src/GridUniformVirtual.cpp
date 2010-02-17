/*
 * GridUniformVirtual.cpp
 *
 *  Created on: Nov 30, 2009
 *      Author: Nico Schl�mer
 */

#include "GridUniformVirtual.h"

// ============================================================================
GridUniformVirtual::GridUniformVirtual( double h,
                                        double gridDomainArea,
                                        unsigned int numGridPoints ) :
GridVirtual( Teuchos::tuple<double>(h,h),
             gridDomainArea,
             numGridPoints )
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
