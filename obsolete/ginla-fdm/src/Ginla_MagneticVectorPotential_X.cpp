#include "Ginla_MagneticVectorPotential_X.h"

// ============================================================================
Ginla::MagneticVectorPotential::X::
X( double mu ) :
  Ginla::MagneticVectorPotential::Virtual( mu )
{
}
// ============================================================================
Ginla::MagneticVectorPotential::X::
~X()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::X::
setParameters( const LOCA::ParameterVector & p )
{
    bool valuesChanged = false;
  
    if (p.isParameter( "H0" ))
        if ( mu_ != p.getValue ( "H0" ) )
        {
            mu_ = p.getValue ( "H0" );
            valuesChanged = true;
        }

    return valuesChanged;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::MagneticVectorPotential::X::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );
          
  p->addParameter( "H0", mu_ );
          
  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::X::
getA(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getAx(x),
                                                          this->getAy(x),
                                                          this->getAz(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::X::
getAx(const Point & x) const
{
  return 0.0;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::X::
getAy(const Point & x) const
{
  return - mu_ * 0.5*x[2];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::X::
getAz(const Point & x) const
{
  return   mu_ * 0.5*x[1];
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::X::
getDADMu(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getDAxDMu(x),
                                                          this->getDAyDMu(x),
                                                          this->getDAzDMu(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::X::
getDAxDMu(const Point & x) const
{
  return 0.0;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::X::
getDAyDMu(const Point & x) const
{
  return - 0.5*x[2];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::X::
getDAzDMu(const Point & x) const
{
  return   0.5*x[1];
}
// ============================================================================