/*
 * MagneticVectorPotential.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: Nico Schlšmer
 */

#include "MagneticVectorPotential.h"

// ============================================================================
MagneticVectorPotential::MagneticVectorPotential( double h0,
                                                  double edgeLength ) :
  h0_(h0),
  edgeLength_(edgeLength)
{
}
// ============================================================================
MagneticVectorPotential::~MagneticVectorPotential()
{
}
// ============================================================================
void
MagneticVectorPotential::setH0( double h0 )
{
  h0_ = h0;
}
// ============================================================================
Teuchos::RCP<Teuchos::Array<double> >
MagneticVectorPotential::getA(const Teuchos::Array<double> & x) const
{
  Teuchos::RCP<Teuchos::Array<double> >  A = Teuchos::rcp( new Teuchos::Array<double>(2) );

  (*A)[0] = getAx(x);
  (*A)[1] = getAy(x);

  return A;
}
// ============================================================================
double
MagneticVectorPotential::getAx(const Teuchos::Array<double> & x) const
{
  return - 0.5 * h0_ * x[1]
         + 0.5 * h0_ * 0.5*edgeLength_;
}
// ============================================================================
double
MagneticVectorPotential::getAy(const Teuchos::Array<double> & x) const
{
  return   0.5 * h0_ * x[0]
         - 0.5 * h0_ * 0.5*edgeLength_;
}
// ============================================================================
