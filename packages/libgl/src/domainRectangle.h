#ifndef DOMAINRECTANGLE_H
#define DOMAINRECTANGLE_H

#include "domainVirtual.h"

class DomainRectangle:
    public DomainVirtual
{
  public:
      //! Default constructor.
      DomainRectangle( double lx, double ly );

      ~DomainRectangle();

      //! This function returns true if the given 2D cartisian coordinate
      //! sits in the domain, and false otherwise.
      bool
      isInDomain( Teuchos::Array<double> & x );

  private:

    double lx_;
    double ly_;

};
#endif // DOMAINRECTANGLE_H