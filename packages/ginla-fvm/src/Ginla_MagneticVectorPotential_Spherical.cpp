#include "Ginla_MagneticVectorPotential_Spherical.h"

// ============================================================================
Ginla::MagneticVectorPotential::Spherical::
Spherical( double mu, double phi, double theta ) :
  Ginla::MagneticVectorPotential::Virtual( mu ),
  phi_( phi ),
  theta_( theta )
{
}
// ============================================================================
Ginla::MagneticVectorPotential::Spherical::
~Spherical()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::Spherical::
setParameters( const LOCA::ParameterVector & p )
{
    bool valuesChanged = false;

    if (p.isParameter( "H0" ))
        if ( mu_ != p.getValue ( "H0" ) )
        {
            mu_ = p.getValue ( "H0" );
            valuesChanged = true;
        }

    if (p.isParameter( "phi" ))
        if ( phi_ != p.getValue ( "phi" ) )
        {
            phi_ = p.getValue ( "phi" );
            valuesChanged = true;
        }

    if (p.isParameter( "theta" ))
        if ( theta_ != p.getValue ( "theta" ) )
        {
            theta_ = p.getValue ( "theta" );
            valuesChanged = true;
        }

    return valuesChanged;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::MagneticVectorPotential::Spherical::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );

  p->addParameter( "H0", mu_ );
  p->addParameter( "phi", phi_ );
  p->addParameter( "theta", theta_ );

  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Spherical::
getA(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getAx(x),
                                                          this->getAy(x),
                                                          this->getAz(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Spherical::
getAx(const Point & x) const
{
  return mu_ * (-0.5 * sin(theta_) * x[1] + 0.5 * cos(theta_) * sin(phi_) * x[2]);
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Spherical::
getAy(const Point & x) const
{
  return  mu_ * (0.5 * sin(theta_) * x[0] - 0.5 * cos(theta_) * cos(phi_) * x[2]);
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Spherical::
getAz(const Point & x) const
{
  return mu_ * ( 0.5 * cos(theta_) * cos(phi_) * x[1] - 0.5 * cos(theta_) * sin(phi_) * x[0] );
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::Spherical::
getDADMu(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getDAxDMu(x),
                                                          this->getDAyDMu(x),
                                                          this->getDAzDMu(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Spherical::
getDAxDMu(const Point & x) const
{
return -0.5 * sin(theta_) * x[1] + 0.5 * cos(theta_) * sin(phi_) * x[2];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Spherical::
getDAyDMu(const Point & x) const
{
  return  0.5 * sin(theta_) * x[0] - 0.5 * cos(theta_) * cos(phi_) * x[2];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::Spherical::
getDAzDMu(const Point & x) const
{
  return 0.5 * cos(theta_) * cos(phi_) * x[1] - 0.5 * cos(theta_) * sin(phi_) * x[0];
}
// ============================================================================