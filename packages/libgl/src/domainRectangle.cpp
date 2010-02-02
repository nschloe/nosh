#include "domainRectangle.h"

// =============================================================================
DomainRectangle::DomainRectangle( double lx, double ly ):
  lx_( lx ),
  ly_( ly )
{
}
// =============================================================================
DomainRectangle::~DomainRectangle()
{
}
// =============================================================================
bool
DomainRectangle::isInDomain( Teuchos::Array<double> & x )
{
    TEST_FOR_EXCEPTION( x.length() != 2,
		                std::logic_error,
		                "Point not 2-dimensional." );

  return ( x[0]>=0.0 && x[0]<=lx_ &&
           x[1]>=0.0 && x[1]<=ly_     );
}
// =============================================================================
