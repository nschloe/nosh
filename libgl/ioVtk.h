#ifndef IOVTK_H
#define IOVTK_H

#include "ioVirtual.h"

#include "staggeredGrid.h"

#include <string>
#include <complex>

#include <Teuchos_ParameterList.hpp>

#include <Tpetra_MultiVector.hpp>

#include <iostream>

typedef std::complex<double> double_complex;

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
    read ( Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > &psi,
           const Teuchos::RCP<const Teuchos::Comm<int> >          comm, // TODO: remove this
           Teuchos::ParameterList                                 &problemParams
         ) const;

    //! Writes the  order parameter \f$\psi\f$ and the problem parameter list
    //! into a legac VTK file.
    //! The data is written such that the lexicographical ordering of the
    //! nodes is preserved and the resulting file contains a state \f$\psi\f
    //! that can be viewed using standard tools.
    virtual void
    write ( const Tpetra::MultiVector<double_complex,int> &psi,
            const Teuchos::ParameterList                  &problemParams,
            const StaggeredGrid                           &sGrid
          ) const;

    virtual void
    write( const Tpetra::MultiVector<double_complex,int> &psi,
           const StaggeredGrid                           &sGrid
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
    writeScalars( const Tpetra::MultiVector<double_complex,int> & psi,
                        const StaggeredGrid                     & sGrid,
                              std::ofstream                     & oStream
                      ) const;

    //! joins a vector of strings to one string with a separator string sep
    std::string strJoin ( const std::vector<std::string> & vec,
                          const std::string              & sep
                        ) const;

  };
#endif // IOVTK_H
