#ifndef IOVIRTUAL_H
#define IOVIRTUAL_H

#include "staggeredGrid.h"

#include <string>
#include <vector>
#include <complex>

#include <Teuchos_ParameterList.hpp>

typedef std::complex<double> double_complex;

class IoVirtual
{
  public:

     //! Default constructor.
     IoVirtual( std::string fname );

     //! Destructor
     virtual ~IoVirtual();

     virtual void read( std::vector<double_complex> *psi,
                        Teuchos::ParameterList      *problemParams ) = 0;

     virtual void write( const std::vector<double_complex> &psi,
                         const Teuchos::ParameterList      &problemParams,
                         StaggeredGrid                     &sGrid          ) = 0;

  protected:
      std::string fileName;

};
#endif // IOVIRTUAL_H