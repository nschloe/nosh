#ifndef IOXDMF_H
#define IOXDMF_H

#include "ioVirtual.h"

#include <Teuchos_XMLObject.hpp>

#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

#include <Tpetra_MultiVector.hpp>

class IoXdmf: public IoVirtual
  {
  public:

    //! Default constructor.
    IoXdmf ( boost::filesystem::path fname );

    //! Destructor
    virtual ~IoXdmf();

    //! Reads the order parameter \f$\psi\f$ and the problem parameter list
    //! from an XDMF file (plus its accompaining HDF5) into the arguments.
    virtual void
    read( const Teuchos::RCP<const Teuchos::Comm<int> > & tComm,
                Teuchos::RCP<DoubleMultiVector>         & x,
                Teuchos::ParameterList                  & problemParams
        ) const;

    //! Virtual function for reading the order parameter \f$\psi\f$ and the
    //! parameter list from a given file.
    virtual void
    read( const Teuchos::RCP<const Teuchos::Comm<int> > & tComm,
                Teuchos::RCP<ComplexMultiVector>        & x,
                Teuchos::ParameterList                  & problemParams
        ) const;

    //! Writes the  order parameter \f$\psi\f$ and the problem parameter list
    //! into an XDMF file (plus its accompaining HDF5).
    //! The data is written such that the lexicographical ordering of the
    //! nodes is preserved and the resulting file contains a state \f$\psi\f
    //! that can be viewed using standard tools.
    virtual void
    write ( const Epetra_MultiVector             & x,
            const Teuchos::Tuple<unsigned int,2> & Nx,
            const Teuchos::Tuple<double,2>       & h,
            const Teuchos::Array<int>            & kBoundingBox,
            const Teuchos::ParameterList         & problemParams,
            const double                         & dummyValue
          );

    virtual void
    write ( const Epetra_MultiVector                        & x,
            const Teuchos::Array<Teuchos::Tuple<double,2> > & loc,
            const Teuchos::ParameterList                    & problemParams
          );
          
    virtual void
    write ( const Epetra_MultiVector             & x,
            const Teuchos::Tuple<unsigned int,2> & Nx,
            const Teuchos::Tuple<double,2>       & h,
            const Teuchos::ParameterList         & problemParams
          );

    virtual void
    write ( const DoubleMultiVector              & x,
            const Teuchos::Tuple<unsigned int,2> & Nx,
            const Teuchos::Tuple<double,2>       & h,
            const Teuchos::ParameterList         & problemParams
          );

    virtual void
    write ( const ComplexMultiVector             & x,
            const Teuchos::Tuple<unsigned int,2> & Nx,
            const Teuchos::Tuple<double,2>       & h,
            const Teuchos::ParameterList         & problemParams
          );

    virtual void
    write ( const ComplexMultiVector              & x,
            const Teuchos::Tuple<unsigned int,2>  & Nx,
            const Teuchos::Tuple<double,2>        & h,
            const Teuchos::Array<int>             & kBoundingBox,
            const Teuchos::ParameterList          & problemParams,
            const double                          & dummyValue
          );


    virtual void
    write( const DoubleMultiVector              & x,
           const Teuchos::Tuple<unsigned int,2> & Nx,
           const Teuchos::Tuple<double,2>       & h
         );

    virtual void
    write( const ComplexMultiVector             & x,
           const Teuchos::Tuple<unsigned int,2> & Nx,
           const Teuchos::Tuple<double,2>       & h
         );

  private:

    //! Read the Heavy Data from a file specified in xmlFileObject.
    //! Note that readVec is given twice-dereferenced. This is because
    //! the crucial part of the routine (the actual HDF5 reader) sets
    //! a pointer to a MultiVector, which we will need to return.
    //! If your pointer is *a, call this routine with &a.
    void getHeavyData ( const Teuchos::XMLObject &xmlFileObject,
                        const Epetra_Comm        &comm,
                        Epetra_MultiVector       **readVec,
                        const std::string        &fileDirectory,
                        const std::string        &xmlName,
                        const std::string        &hdf5GroupName
                      ) const;

    const Teuchos::XMLObject* xmlFind ( const Teuchos::XMLObject *xmlObj,
                                        const std::string        tag
                                      ) const;

    const Teuchos::XMLObject* xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
                                                 const std::string        tag,
                                                 const std::string        attribute,
                                                 const std::string        value
                                               ) const;

  };
#endif // IOXDMF_H
