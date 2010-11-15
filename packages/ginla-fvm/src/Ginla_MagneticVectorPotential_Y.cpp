#include "Ginla_MagneticVectorPotential_Y.h"

// ============================================================================
Ginla::MagneticVectorPotential::Y::
Y( double mu ) :
  mu_( mu )
{
}
// ============================================================================
Ginla::MagneticVectorPotential::Y::
~Y()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::Y::
setParameters( const LOCA::ParameterVector & p )
{
    bool valuesChanged = false;

    if (p.isParameter( "mu" ))
        if ( mu_ != p.getValue ( "mu" ) )
        {
            mu_ = p.getValue ( "mu" );
            valuesChanged = true;
        }

    return valuesChanged;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::MagneticVectorPotential::Y::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );

  p->addParameter( "mu", mu_ );

  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Y::
getA(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getAx(x),
                                                          this->getAy(x),
                                                          this->getAz(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Y::
getAx(const Point & x) const
{
  return   mu_ * 0.5*x[2];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Y::
getAy(const Point & x) const
{
  return 0.0;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Y::
getAz(const Point & x) const
{
  return - mu_ * 0.5*x[0];
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Y::
getDADMu(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getDAxDMu(x),
                                                          this->getDAyDMu(x),
                                                          this->getDAzDMu(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Y::
getDAxDMu(const Point & x) const
{
  return 0.5*x[2];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Y::
getDAyDMu(const Point & x) const
{
  return 0.0;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Y::
getDAzDMu(const Point & x) const
{
  return - 0.5*x[0];
}
// ============================================================================