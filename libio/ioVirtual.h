#ifndef IOVIRTUAL_H
#define IOVIRTUAL_H

#include <string>

#include <Teuchos_ParameterList.hpp>

#include <Tpetra_MultiVector.hpp>

#include <Teuchos_Comm.hpp>
#include <Teuchos_Tuple.hpp>

#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

typedef Tpetra::MultiVector<double,Thyra::Ordinal> DoubleMultiVector;

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
    read( const Teuchos::RCP<const Teuchos::Comm<int> > & tComm,
                Teuchos::RCP<DoubleMultiVector>         & x,
                Teuchos::ParameterList                  & problemParams
        ) const = 0; // pure virtual

    //! Virtual function for writing the order parameter \f$\psi\f$ and the
    //! parameter list to a given file.
    virtual void
    write ( const DoubleMultiVector              & x,
            const Teuchos::Tuple<unsigned int,2> & Nx,
            const Teuchos::Tuple<double,2>       & h,
            const Teuchos::ParameterList         & problemParams
          ) = 0; // pure virtual

    virtual void
    write( const DoubleMultiVector              & x,
           const Teuchos::Tuple<unsigned int,2> & Nx,
           const Teuchos::Tuple<double,2>       & h
          ) = 0; // pure virtual

  protected:
    //! File name for the I/O.
    std::string fileName_;

  };
#endif // IOVIRTUAL_H
