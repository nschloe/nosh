#ifndef IOVIRTUAL_H
#define IOVIRTUAL_H



#include "staggeredGrid.h"

#include <string>
#include <complex>

#include <Teuchos_ParameterList.hpp>

#include <Tpetra_MultiVector.hpp>

typedef std::complex<double> double_complex;

class IoVirtual
  {
  public:

    //! Default constructor.
    IoVirtual ( std::string fname );

    //! Destructor
    virtual ~IoVirtual();

    //! Virtual function for reading the order parameter \f$\psi\f$ and the
    //! parameter list from a given file.
    virtual void
    read ( Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi,
           Teuchos::RCP<Teuchos::Comm<int> >                      comm, // TODO: remove this
           Teuchos::ParameterList                                 *problemParams
         ) = 0; // pure virtual

    //! Virtual function for writing the order parameter \f$\psi\f$ and the
    //! parameter list to a given file.
    virtual void
    write ( const Tpetra::MultiVector<double_complex,int> &psi,
            const Teuchos::ParameterList                  &problemParams,
            StaggeredGrid                                 &sGrid
          ) = 0; // pure virtual

  protected:
    //! File name for the I/O.
    std::string fileName;

  };
#endif // IOVIRTUAL_H
