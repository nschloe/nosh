/*
 * MagneticVectorPotential.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: Nico Schl\"omer
 */

#include "Ginla_MagneticVectorPotential_Centered.h"

// ============================================================================
Ginla::MagneticVectorPotential::Centered::
Centered( double h0,
          double edgeLength ) :
  h0_(h0),
  edgeLength_(edgeLength)
{
}
// ============================================================================
Ginla::MagneticVectorPotential::Centered::
~Centered()
{
}
// ============================================================================
void
Ginla::MagneticVectorPotential::Centered::
setParameters( const LOCA::ParameterVector & p)
{
    h0_ = p.getValue ( "H0" );
    edgeLength_ = p.getValue ( "scaling" );
    return;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::MagneticVectorPotential::Centered::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );
          
  p->addParameter( "H0", h0_ );
  p->addParameter( "scaling", edgeLength_ );
          
  return p;
}
// ============================================================================
Teuchos::RCP<Teuchos::Array<double> >
Ginla::MagneticVectorPotential::Centered::
getA(const Teuchos::Array<double> & x) const
{
  Teuchos::RCP<Teuchos::Array<double> >  A = Teuchos::rcp( new Teuchos::Array<double>(2) );

  (*A)[0] = getAx(x);
  (*A)[1] = getAy(x);

  return A;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getAx(const Teuchos::Array<double> & x) const
{
  return - 0.5 * h0_ * x[1]
         + 0.5 * h0_ * 0.5*edgeLength_;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getAy(const Teuchos::Array<double> & x) const
{
  return   0.5 * h0_ * x[0]
         - 0.5 * h0_ * 0.5*edgeLength_;
}
// ============================================================================