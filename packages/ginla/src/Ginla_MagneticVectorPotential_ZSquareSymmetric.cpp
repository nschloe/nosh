#include "Ginla_MagneticVectorPotential_ZSquareSymmetric.h"



// ============================================================================
Ginla::MagneticVectorPotential::ZSquareSymmetric::
ZSquareSymmetric( double mu,
                  double edgeLength ) :
  Ginla::MagneticVectorPotential::Virtual( mu ),
  edgeLength_( edgeLength )
{
}
// ============================================================================
Ginla::MagneticVectorPotential::ZSquareSymmetric::
~ZSquareSymmetric()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::ZSquareSymmetric::
setParameters( const LOCA::ParameterVector & p )
{
    bool valuesChanged = false;
  
    if (p.isParameter( "H0" ))
        if ( mu_ != p.getValue ( "H0" ) )
        {
            mu_ = p.getValue ( "H0" );
            valuesChanged = true;
        }
        
    if (p.isParameter( "edge length" ))
        if ( edgeLength_ != p.getValue ( "edge length" ) )
        {
            edgeLength_ = p.getValue ( "edge length" );
            valuesChanged = true;
        }

    return valuesChanged;
}
// ============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::MagneticVectorPotential::ZSquareSymmetric::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );
          
  p->addParameter( "H0", mu_ );
  p->addParameter( "edge length", edgeLength_ );
          
  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::ZSquareSymmetric::
getA(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getAx(x),
                                                          this->getAy(x),
                                                          this->getAz(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::ZSquareSymmetric::
getAx(const Point & x) const
{
  return - 0.5 * mu_ * x[1]
         + 0.5 * mu_ * 0.5*edgeLength_;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::ZSquareSymmetric::
getAy(const Point & x) const
{
  return   0.5 * mu_ * x[0]
         - 0.5 * mu_ * 0.5*edgeLength_;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::ZSquareSymmetric::
getAz(const Point & x) const
{
  return 0.0;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::ZSquareSymmetric::
getDADMu(const Point & x) const
{
  return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getDAxDMu(x),
                                                          this->getDAyDMu(x),
                                                          this->getDAzDMu(x)
                     ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::ZSquareSymmetric::
getDAxDMu(const Point & x) const
{
  return - 0.5 * x[1]
         + 0.5 * 0.5*edgeLength_;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::ZSquareSymmetric::
getDAyDMu(const Point & x) const
{
  return   0.5 * x[0]
         - 0.5 * 0.5*edgeLength_;
}
// ============================================================================
double
Ginla::MagneticVectorPotential::ZSquareSymmetric::
getDAzDMu(const Point & x) const
{
  return 0.0;
}
// ============================================================================