#include "stateFileReader.h"

#include <Teuchos_FileInputSource.hpp>
#include <Teuchos_XMLParameterListReader.hpp>

// =============================================================================
// Constructor
StateFileReader::StateFileReader()
{
}
// =============================================================================



// =============================================================================
// Destructor
StateFileReader::~StateFileReader()
{
}
// =============================================================================



// =============================================================================
void StateFileReader::readFile( std::string                 &fileName,
                                std::string                 &fileFormat,
                                std::vector<double_complex> *psi,
                                Teuchos::ParameterList      *problemParams )
{
  if ( !fileFormat.compare("legacyVTK") )
      readLegacyVtkFile( fileName, psi, problemParams );
  else if ( !fileFormat.compare("VTI")
 )
      readVtiFile( fileName, psi, problemParams );
  else if ( !fileFormat.compare("XDMF") )
      readXdmfFile( fileName, psi, problemParams );
  else {
      std::cerr << "StateFileReader::readFile :" << std::endl
                << "Illegal file format \"" << fileFormat <<"\"." << std::endl
                << "Use one of \"legacyVTK\", \"VTI\", and \"XDMF\"." << std::endl;
      exit(EXIT_FAILURE);
  }
}
// =============================================================================



// =============================================================================
void StateFileReader::readLegacyVtkFile( std::string                 &filename,
                                         std::vector<double_complex> *psi,
                                         Teuchos::ParameterList      *problemParams )
{
  std::cerr << "StateFileReader::readLegacyVtkFile not yet implemented." << std::endl;
  exit(EXIT_FAILURE);
}
// =============================================================================




// =============================================================================
void StateFileReader::readVtiFile( std::string                 &filename,
                                   std::vector<double_complex> *psi,
                                   Teuchos::ParameterList      *problemParams )
{

  std::cerr << "StateFileReader::readVtiFile not yet implemented." << std::endl;
  exit(EXIT_FAILURE);

  // pass a possible 
  //<?xml version="1.0"?>
  // at the beginning of the file

  std::cout << "11" << std::endl;
    Teuchos::FileInputSource xmlFile(filename);

  std::cout << filename << std::endl;

    // Extract the object from the filename.
    // -- This is actually quite costly, as it read all -- *all* -- data in as
    // strings, including the heavy numerical data.
    Teuchos::XMLObject xmlFileObject = xmlFile.getObject();

    // find and read the parameter list
    const Teuchos::XMLObject* parameterListObject =
                                      xmlFind( &xmlFileObject, "ParameterList" );
//     problemParams = Teuchos::XMLParameterListReader().
//                                           toParameterList( *parameterListObject );

    // find and get the absolute value vector
    std::string absPsiString;
    const Teuchos::XMLObject* absPsiObject = xmlAttributeFind( &xmlFileObject,
                                                              "DataArray",
                                                              "Name",
                                                              "abs(psi)"      );

    if ( !absPsiObject ) { // pointer invalid
        std::cerr << "No such XML Object found! Abort." << std::endl;
        exit(EXIT_FAILURE);
    }
  
  std::cout << "1a" << std::endl;
  std::cout << *absPsiObject << std::endl;
  std::cout << "1b" << std::endl;
  std::cout << absPsiObject->numContentLines() << std::endl;
  std::cout << "1c" << std::endl;
  
  
  //   const char* str;
  //   char** endptr;
  //   for (int l=0; l<absPsiObject->numContentLines(); l++) {
  //       // get the full row into a char*
  //       str = absPsiObject->getContentLine(l).c_str();
  // std::cout << "l = " << l << std::endl;
  // std::cout << "line: " << absPsiObject->getContentLine(l) << std::endl;
  // std::cout << str << std::endl;
  // std::cout << strtod( str , endptr ) << std::endl;
  //   }
  
  // std::cout << "absPsiString :" << absPsiString << std::endl;
  
  //   // plot the contents
  //   std::cout << xmlFileObject << std::endl;
  
  std::cout << "22" << std::endl;
  
    std::cout << " plist: " << problemParams << std::endl;
  
  std::cout << "33" << std::endl;

}
// =============================================================================




// =============================================================================
void StateFileReader::readXdmfFile( std::string                 &filename,
                                    std::vector<double_complex> *psi,
                                    Teuchos::ParameterList      *problemParams )
{

  std::cerr << "StateFileReader::readXdmfFile not yet implemented." << std::endl;
  exit(EXIT_FAILURE);

}
// =============================================================================




// =============================================================================
// Inside an XML object, this function looks for a specific tag and returns
// a pointer to it.
const Teuchos::XMLObject* StateFileReader::xmlFind ( const Teuchos::XMLObject *xmlObj,
                                                     const std::string        tag )
{
  const Teuchos::XMLObject* xmlOut=NULL;

  if ( !xmlObj->getTag().compare(tag) ) // strings are equal
      xmlOut = xmlObj;
  else
      for (int k=0; k<xmlObj->numChildren(); k++) {
          xmlOut = xmlFind ( &(xmlObj->getChild(k)), tag ); // recursive call
          if (xmlOut) break; // not the null pointer => return
      }

  return xmlOut;
}
// =============================================================================



// =============================================================================
const Teuchos::XMLObject* StateFileReader::xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
                                                              const std::string        tag,
                                                              const std::string        attribute,
                                                              const std::string        value      )
{
  const Teuchos::XMLObject* xmlOut=NULL;

// std::cout << "Looking for: " << tag << " " << attribute << " " <<  value << std::endl;
// std::cout << "Found:       " << xmlObj->getTag() << " " << xmlObj->hasAttribute(attribute) << std::endl;
// if ( xmlObj->hasAttribute(attribute) )
//      std::cout << "1" << std::endl;
//      std::cout << xmlObj->getAttribute(attribute) << std::endl;

  if (    !xmlObj->getTag().compare(tag)
       && xmlObj->hasAttribute(attribute)
       && xmlObj->getAttribute(attribute).compare(value) )
      xmlOut = xmlObj;
  else
      for (int k=0; k<xmlObj->numChildren(); k++) {
          xmlOut = xmlAttributeFind ( &(xmlObj->getChild(k)), // recursive call
                                      tag,
                                      attribute,
                                      value                    );
          if (xmlOut) break; // not the null pointer => return
      }

  return xmlOut;
}
// =============================================================================