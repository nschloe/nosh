#include "domainRectangle.h"

#include "glException.h"

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
  if ( x.length() != 2 ) {
      throw glException( "DomainRectangle::isInDomain",
                         "Point not 2-dimensional." );
  }

  return ( x[0]>=0.0 && x[0]<=lx_ &&
           x[1]>=0.0 && x[1]<=ly_     );
}
// =============================================================================