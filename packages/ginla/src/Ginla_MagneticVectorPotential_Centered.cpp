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
bool
Ginla::MagneticVectorPotential::Centered::
setParameters( const LOCA::ParameterVector & p )
{
    bool valuesChanged = false;
  
    if (p.isParameter( "H0" ))
        if ( h0_ != p.getValue ( "H0" ) )
        {
            h0_ = p.getValue ( "H0" );
            valuesChanged = true;
        }
           
    
    if ( p.isParameter("scaling") )
        if ( edgeLength_ != p.getValue ( "scaling" ) )
        {
            edgeLength_ = p.getValue ( "scaling" );
            valuesChanged = true;
        }

    return valuesChanged;
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
Teuchos::RCP<DoubleTuple>
Ginla::MagneticVectorPotential::Centered::
getA(const DoubleTuple & x) const
{
  return Teuchos::rcp( new DoubleTuple( Teuchos::tuple<double>( this->getAx(x),
                                                                this->getAy(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getAx(const DoubleTuple & x) const
{
  return - 0.5 * h0_ * x[1]
         + 0.5 * h0_ * 0.5*edgeLength_;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getAy(const DoubleTuple & x) const
{
  return   0.5 * h0_ * x[0]
         - 0.5 * h0_ * 0.5*edgeLength_;
}
// ============================================================================
Teuchos::RCP<DoubleTuple>
Ginla::MagneticVectorPotential::Centered::
getDADh0(const DoubleTuple & x) const
{
  return Teuchos::rcp( new DoubleTuple( Teuchos::tuple<double>( this->getDAxDh0(x),
                                                                this->getDAyDh0(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getDAxDh0(const DoubleTuple & x) const
{
  return - 0.5 * x[1]
         + 0.5 * 0.5*edgeLength_;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getDAyDh0(const DoubleTuple & x) const
{
  return   0.5 * x[0]
         - 0.5 * 0.5*edgeLength_;
}
// ============================================================================