#include "Ginla_MagneticVectorPotential_MagneticDot.h"

// ============================================================================
Ginla::MagneticVectorPotential::MagneticDot::
MagneticDot( double mu ) :
  Ginla::MagneticVectorPotential::Virtual( mu ),
  magnetRadius_( 4.5315 ),
  zz1_ ( 0.1 ),
  zz2_ ( 1.1 )
{
}
// =============================================================================
Ginla::MagneticVectorPotential::MagneticDot::
~MagneticDot()
{
}
// ============================================================================
bool
Ginla::MagneticVectorPotential::MagneticDot::
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
Ginla::MagneticVectorPotential::MagneticDot::
getParameters() const
{
  Teuchos::RCP<LOCA::ParameterVector> p =
          Teuchos::rcp( new LOCA::ParameterVector() );

  p->addParameter( "mu", mu_ );

  return p;
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::MagneticDot::
getA( const Point & x ) const
{
    // Span a cartesian grid over the sample, and integrate over it.

    // The number of grids within the diameter of the magnetic dot, counting
    // the boundary nodes in.
    int n_mag = 101;
    // grid width
    double dx = 2.0 * magnetRadius_ / (n_mag - 1);

    double ax = 0.0;
    double ay = 0.0;

    for ( int ix=0; ix < n_mag; ix++ )
    {
        double xi = ix*dx - magnetRadius_ ;
        for ( int iy=0; iy < n_mag; iy++ )
        {
            double yi = iy*dx - magnetRadius_;
            // circular shape magnetic dot
            if ( xi*xi + yi*yi <= magnetRadius_*magnetRadius_ )
            {
                // x distance between grid point x to magnetic point xi
                double xx = x[0] - xi;
                // y distance between grid point y to magnetic point yi
                double yy = x[1] - yi;
                // r distance between grid point X to magnetic point (xi,yi)
                double r = xx * xx + yy * yy;

                if ( r > 1.0e-15 )
                {
                    // 3D distance to point on upper edge (xi,yi,zz1)
                    double r_3D1 = sqrt( r + zz1_ * zz1_);
                    // 3D distance to point on lower edge (xi,yi,zz2)
                    double r_3D2 = sqrt( r + zz2_ * zz2_);
                    double alpha = ( zz2_ / r_3D2 - zz1_ / r_3D1 ) * dx * dx / r;
                    ax += yy * alpha;
                    ay -= xx * alpha;
                }
            }
        }
    }

    ax *= mu_;
    ay *= mu_;

    return Teuchos::rcp( new Point( Teuchos::tuple<double>( ax, ay, 0.0 ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::MagneticDot::
getAx(const Point & x) const
{
    Teuchos::RCP<Point> a = this->getA( x );
    return (*a)[0];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::MagneticDot::
getAy(const Point & x) const
{
    Teuchos::RCP<Point> a = this->getA( x );
    return (*a)[1];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::MagneticDot::
getAz(const Point & x) const
{
    Teuchos::RCP<Point> a = this->getA( x );
    return (*a)[2];
}
// ============================================================================
Teuchos::RCP<Point>
Ginla::MagneticVectorPotential::MagneticDot::
getDADMu(const Point & x) const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Not yet implemented." );
    return Teuchos::rcp( new Point( Teuchos::tuple<double>( this->getDAxDMu(x),
                                                            this->getDAyDMu(x),
                                                            this->getDAzDMu(x)
                                                            ) ) );
}
// ============================================================================
double
Ginla::MagneticVectorPotential::MagneticDot::
getDAxDMu(const Point & x) const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Not yet implemented." );
    return - 0.5 * x[1];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::MagneticDot::
getDAyDMu(const Point & x) const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Not yet implemented." );
    return 0.5 * x[0];
}
// ============================================================================
double
Ginla::MagneticVectorPotential::MagneticDot::
getDAzDMu(const Point & x) const
{
    return 0.0;
}
// ============================================================================
