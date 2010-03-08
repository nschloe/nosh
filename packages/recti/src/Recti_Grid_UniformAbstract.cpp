/*
 * GridUniformVirtual.cpp
 *
 *  Created on: Nov 30, 2009
 *      Author: Nico Schl\"omer
 */

#include "Recti_Grid_UniformAbstract.h"

// ============================================================================
Recti::Grid::UniformAbstract::
UniformAbstract( double h,
                    double gridDomainArea,
                    unsigned int numGridPoints ) :
Abstract( Teuchos::tuple<double>(h,h),
             gridDomainArea,
             numGridPoints )
{
}
// ============================================================================
Recti::Grid::UniformAbstract::
UniformAbstract() :
Abstract()
{
}
// ============================================================================
Recti::Grid::UniformAbstract::
~UniformAbstract()
{
}
// ============================================================================
double
Recti::Grid::UniformAbstract::
getUniformH() const
{
  return h_[0];
}
// ============================================================================
