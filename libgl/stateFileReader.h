#include <string>
#include <vector>
#include <complex>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLObject.hpp>

typedef std::complex<double> double_complex;

class StateFileReader
{
  public:

     //! Default constructor. 
     StateFileReader();

     //! Destructor
     ~StateFileReader();

     //
     void readFile( std::string                 &fileName,
                    std::string                 &fileFormat,
                    std::vector<double_complex> *psi,
                    Teuchos::ParameterList      *problemParams );

  private:

     void readLegacyVtkFile( std::string                 &filename,
                             std::vector<double_complex> *psi,
                             Teuchos::ParameterList      *problemParams );
     void readVtiFile      ( std::string                 &filename,
                             std::vector<double_complex> *psi,
                             Teuchos::ParameterList      *problemParams );
     void readXdmfFile     ( std::string                 &filename,
                             std::vector<double_complex> *psi,
                             Teuchos::ParameterList      *problemParams );

     const Teuchos::XMLObject* xmlFind ( const Teuchos::XMLObject *xmlObj,
                                         const std::string        tag );

     const Teuchos::XMLObject* xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
                                                  const std::string        tag,
                                                  const std::string        attribute,
                                                  const std::string        value      );
};