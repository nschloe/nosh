#ifndef DOMAINVIRTUAL_H
#define DOMAINVIRTUAL_H

#include <Teuchos_Array.hpp>

class DomainVirtual
{
  public:
      //! Default constructor.
      DomainVirtual();

      virtual ~DomainVirtual();

      //! This function returns true if the given 2D cartisian coordinate
      //! sits in the domain, and false otherwise.
      virtual bool
      isInDomain( Teuchos::Array<double> & x ) = 0;

};
#endif // DOMAINVIRTUAL_H