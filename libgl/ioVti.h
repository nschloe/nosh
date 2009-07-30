#include "staggeredGrid.h"

#include <string>
#include <vector>
#include <complex>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLObject.hpp>

typedef std::complex<double> double_complex;

class IoVti
{
  public:

     //! Default constructor.
     IoVti( StaggeredGrid &sGrid );

     //! Destructor
     ~IoVti();

     void read( const std::string           &fileName,
                std::vector<double_complex> *psi,
                Teuchos::ParameterList      *problemParams );

     //
     void write( const std::string                 &fileName,
                 const std::vector<double_complex> &psi,
                 const Teuchos::ParameterList      &problemParams );

  private:

     StaggeredGrid::StaggeredGrid SGrid;

     const Teuchos::XMLObject* xmlFind ( const Teuchos::XMLObject *xmlObj,
                                         const std::string        tag );

     const Teuchos::XMLObject* xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
                                                  const std::string        tag,
                                                  const std::string        attribute,
                                                  const std::string        value      );

};