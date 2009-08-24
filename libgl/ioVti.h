#include "ioVirtual.h"

#include <Teuchos_XMLObject.hpp>

class IoVti: public IoVirtual
{
  public:

     //! Default constructor.
     IoVti( std::string fname );

     //! Destructor
     virtual ~IoVti();

     //! Reads the order parameter \f$\psi\f$ and the problem parameter list
     //! from a VTI file into the arguments.
     virtual void read( std::vector<double_complex> *psi,
                        Teuchos::ParameterList      *problemParams );

     //! Writes the  order parameter \f$\psi\f$ and the problem parameter list
     //! into an XML-style VTI file.
     //! The data is written such that the lexicographical ordering of the
     //! nodes is preserved and the resulting file contains a state \f$\psi\f
     //! that can be viewed using standard tools.
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