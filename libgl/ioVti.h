#include "ioVirtual.h"

#include <Teuchos_XMLObject.hpp>

class IoVti: public IoVirtual
{
  public:

     //! Default constructor.
     IoVti( std::string fname );

     //! Destructor
     virtual ~IoVti();

     virtual void read( std::vector<double_complex> *psi,
                        Teuchos::ParameterList      *problemParams );

     virtual void write( const std::vector<double_complex> &psi,
                         const Teuchos::ParameterList      &problemParams,
                         StaggeredGrid                     &sGrid          );

  private:

     const Teuchos::XMLObject* xmlFind ( const Teuchos::XMLObject *xmlObj,
                                         const std::string        tag );

     const Teuchos::XMLObject* xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
                                                  const std::string        tag,
                                                  const std::string        attribute,
                                                  const std::string        value      );

};