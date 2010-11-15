#include "Ginla_MagneticVectorPotential_Z.h"

// ============================================================================
Ginla::MagneticVectorPotential::Z::
Z( double mu ) :
  mu_( mu )
{
}
// ============================================================================
Ginla::MagneticVectorPotential::Z::
~Z()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::Z::
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
Ginla::MagneticVectorPotential::Z::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );

  p->addParameter( "mu", mu_ );

  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Z::
getA(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getAx(x),
                                                          this->getAy(x),
                                                          this->getAz(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Z::
getAx(const Point & x) const
{
  return - 0.5 * mu_ * x[1];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Z::
getAy(const Point & x) const
{
  return   0.5 * mu_ * x[0];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Z::
getAz(const Point & x) const
{
  return 0.0;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Z::
getDADMu(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getDAxDMu(x),
                                                          this->getDAyDMu(x),
                                                          this->getDAzDMu(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Z::
getDAxDMu(const Point & x) const
{
  return - 0.5 * x[1];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Z::
getDAyDMu(const Point & x) const
{
  return 0.5 * x[0];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Z::
getDAzDMu(const Point & x) const
{
  return 0.0;
}
// ============================================================================
