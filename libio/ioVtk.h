#ifndef IOVTK_H
#define IOVTK_H

#include "ioVirtual.h"

#include <string>
#include <complex>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_MultiVector.hpp>

#include <iostream>

class IoVtk: public IoVirtual
  {
  public:

    //! Default constructor.
    IoVtk ( std::string fname );

    //! Destructor
    virtual ~IoVtk();

    //! Reads the order parameter \f$\psi\f$ and the problem parameter list
    //! from a legacy VTK file into the arguments.
    virtual void
    read( const Teuchos::RCP<const Teuchos::Comm<int> >        &tComm,
                Teuchos::RCP<Tpetra::MultiVector<double,int> > &x,
                Teuchos::ParameterList                         &problemParams
        ) const;

    //! Writes the  order parameter \f$\psi\f$ and the problem parameter list
    //! into a legac VTK file.
    //! The data is written such that the lexicographical ordering of the
    //! nodes is preserved and the resulting file contains a state \f$\psi\f
    //! that can be viewed using standard tools.
    virtual void
    write ( const Tpetra::MultiVector<double,int> & x,
            const int                               Nx,
            const double                            h,
            const Teuchos::ParameterList          & problemParams
          ) const;

    virtual void
    write( const Tpetra::MultiVector<double,int> & x,
           const int                               Nx,
           const double                            h
         ) const;

  protected:
  private:

    void
    writeVtkStructuredPointsHeader( std::ofstream & ioStream,
                                   const int     Nx,
                                   const double  h,
                                   const int     numScalars
                                 ) const;
    
    void
    writeParameterList( const Teuchos::ParameterList & pList,
                        std::ofstream                & ioStream
                      ) const;

    void
    writeScalars( const Tpetra::MultiVector<double,int> & psi,
                  const int                               Nx,
                        std::ofstream                   & oStream
                ) const;

    //! joins a vector of strings to one string with a separator string sep
    std::string strJoin ( const std::vector<std::string> & vec,
                          const std::string              & sep
                        ) const;
  
    bool
    ReadParamsFromVtkFile( std::ifstream          &iFile,
                           Teuchos::ParameterList &fileParams
                         ) const;

    int
    readVtkHeader( std::ifstream & iFile
                 ) const;

    bool
    ReadScalarsFromVtkFile( std::ifstream                                  & iFile,
                            Teuchos::RCP<Tpetra::MultiVector<double,int> > & scalars
                          ) const;

  };
#endif // IOVTK_H
