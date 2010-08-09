#include "Ginla_MagneticVectorPotential_Centered.h"



// ============================================================================
Ginla::MagneticVectorPotential::Centered::
Centered( double h0 ) :
  h0_(h0)
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
setH0( const double mu )
{
  bool valuesChanged = false;

  if ( h0_ != mu )
  {
      h0_ = mu;
      valuesChanged = true;
  }

  return valuesChanged;
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
          
  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Centered::
getA(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getAx(x),
                                                          this->getAy(x),
                                                          this->getAz(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getAx(const Point & x) const
{
  return - 0.5 * h0_ * x[1];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getAy(const Point & x) const
{
  return   0.5 * h0_ * x[0];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getAz(const Point & x) const
{
  return 0.0;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Centered::
getDADh0(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getDAxDh0(x),
                                                          this->getDAyDh0(x),
                                                          this->getDAzDh0(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getDAxDh0(const Point & x) const
{
  return - 0.5 * x[1];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getDAyDh0(const Point & x) const
{
  return 0.5 * x[0];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Centered::
getDAzDh0(const Point & x) const
{
  return 0.0;
}
// ============================================================================