#ifndef IOVTK_H
#define IOVTK_H

#include "ioVirtual.h"

#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_MultiVector.hpp>

#include <iostream>

class IoVtk : public IoVirtual
{
public:

    //! Default constructor.
    IoVtk ( std::string fname );

    //! Destructor
    virtual
    ~IoVtk();

    //! Reads the order parameter \f$\psi\f$ and the problem parameter list
    //! from a legacy VTK file into the arguments.
    virtual void
    read ( const Teuchos::RCP<const Teuchos::Comm<int> > &tComm,
           Teuchos::RCP<DoubleMultiVector> &x,
           Teuchos::ParameterList &problemParams ) const;

    //! Reads the order parameter \f$\psi\f$ and the problem parameter list
    //! from a legacy VTK file into the arguments.
    virtual void
    read ( const Teuchos::RCP<const Teuchos::Comm<int> > &tComm,
           Teuchos::RCP<ComplexMultiVector> &x,
           Teuchos::ParameterList &problemParams ) const;

    //! Writes the  order parameter \f$\psi\f$ and the problem parameter list
    //! into a legacy VTK file.
    //! The data is written such that the lexicographical ordering of the
    //! nodes is preserved and the resulting file contains a state \f$\psi\f
    //! that can be viewed using standard tools.
    virtual void
    write ( const Epetra_MultiVector             & x,
            const Teuchos::Tuple<unsigned int,2> & Nx,
            const Teuchos::Tuple<double,2>       & h,
            const Teuchos::Array<int>            & filter,
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
    write ( const DoubleMultiVector               & x,
            const Teuchos::Tuple<unsigned int,2>  & Nx,
            const Teuchos::Tuple<double,2>        & h,
            const Teuchos::ParameterList          & problemParams );

    virtual void
    write ( const ComplexMultiVector              & x,
            const Teuchos::Tuple<unsigned int,2>  & Nx,
            const Teuchos::Tuple<double,2>        & h,
            const Teuchos::ParameterList          & problemParams );

    virtual void
    write ( const DoubleMultiVector               & x,
            const Teuchos::Tuple<unsigned int,2>  & Nx,
            const Teuchos::Tuple<double,2>        & h
          );

    virtual void
    write ( const ComplexMultiVector             & x,
            const Teuchos::Tuple<unsigned int,2> & Nx,
            const Teuchos::Tuple<double,2>       & h
          );

protected:
private:

    void
    writeVtkStructuredPointsHeader ( std::ofstream & ioStream,
                                     const Teuchos::Tuple<unsigned int,2>  & Nx,
                                     const Teuchos::Tuple<double,2>        & h,
                                     const int numScalars
                                   ) const;

    void
    writeVtkStructuredGridHeader ( std::ofstream                         & ioStream,
                                   const Teuchos::Tuple<unsigned int,2>  & Nx
                                 ) const;

    void
    writeParameterList ( const Teuchos::ParameterList & pList,
                         std::ofstream & ioStream ) const;


    void
    writePointsData ( std::ofstream                                   & ioStream,
                      const Teuchos::Array<Teuchos::Tuple<double,2> > & loc
                    ) const;

    void
    writeScalarsPointData ( const Epetra_MultiVector & x,
                            std::ofstream            & oStream
                          ) const;

    void
    writeScalars ( const Epetra_MultiVector  & x,
                   const Teuchos::Array<int> & filter,
                   const double                dummyValue,
                   std::ofstream             & oStream
                 ) const;

    void
    writeScalars ( const Epetra_MultiVector & x,
                   std::ofstream            & oStream
                 ) const;

    void
    writeScalars ( const DoubleMultiVector  & x,
                   std::ofstream            & oStream
                 ) const;
    void
    writeScalars ( const ComplexMultiVector & psi,
                   std::ofstream            & oStream
                 ) const;

    //! joins a vector of strings to one string with a separator string \c sep
    std::string
    strJoin ( const std::vector<std::string> & vec, const std::string & sep ) const;

    void
    ReadParamsFromVtkFile ( const std::string            & aString,
                            Teuchos::ParameterList & fileParams
                          ) const;

    void
    readVtkHeader ( std::ifstream          & iFile,
                    int                    * vecSize,
                    Teuchos::ParameterList * paramList
                  ) const;

    void
    ReadScalarsFromVtkFile ( std::ifstream                   & iFile,
                             const unsigned int                pointData,
                             Teuchos::RCP<DoubleMultiVector> & scalars
                           ) const;

    void
    ReadScalarsFromVtkFile ( std::ifstream                    & iFile,
                             const unsigned int                 pointData,
                             Teuchos::RCP<ComplexMultiVector> & scalars
                           ) const;

    void
    readScalarFieldHeader ( std::ifstream & iFile,
                            std::string   & buf,
                            int           & numComponents
                          ) const;

};
#endif // IOVTK_H
