#include "ioVirtual.h"

#include <Teuchos_XMLObject.hpp>

#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

class IoXdmf: public IoVirtual
{
  public:

     //! Default constructor.
     IoXdmf( std::string fname );

     //! Destructor
     virtual ~IoXdmf();

     virtual void read( std::vector<double_complex> *psi,
                        Teuchos::ParameterList      *problemParams );

     virtual void write( const std::vector<double_complex> &psi,
                         const Teuchos::ParameterList      &problemParams,
                         StaggeredGrid                     &sGrid          ) ;

  private:

      //! Read the Heavy Data from a file specified in xmlFileObject.
      //! Note that readVec is given twice-dereferenced. This is because
      //! the crucial part of the routine (the actual HDF5 reader) sets
      //! a pointer to a MultiVector, which we will need to return.
      //! If your pointer is *a, call this routine with &a.
      void getHeavyData( const Teuchos::XMLObject &xmlFileObject,
                         const Epetra_Comm        &comm,
                         Epetra_MultiVector       **readVec,
                         const std::string        &fileDirectory,
                         const std::string        &xmlName,
                         const std::string        &hdf5GroupName );

      inline double_complex polar2complex( double abs,
                                           double arg  );

      const Teuchos::XMLObject* xmlFind ( const Teuchos::XMLObject *xmlObj,
                                          const std::string        tag      );

      const Teuchos::XMLObject* xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
                                                   const std::string        tag,
                                                   const std::string        attribute,
                                                   const std::string        value      );

};