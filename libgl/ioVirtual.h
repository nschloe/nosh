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

     //! Virtual function for reading the order parameter \f$\psi\f$ and the
     //! parameter list from a given file.
     virtual void read( std::vector<double_complex> *psi,
                        Teuchos::ParameterList      *problemParams ) = 0;

     //! Virtual function for writing the order parameter \f$\psi\f$ and the
     //! parameter list to a given file.
     virtual void write( const std::vector<double_complex> &psi,
                         const Teuchos::ParameterList      &problemParams,
                         StaggeredGrid                     &sGrid          ) = 0;

  protected:
      //! File name for the I/O.
      std::string fileName;

};
#endif // IOVIRTUAL_H