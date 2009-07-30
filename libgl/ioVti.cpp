#include "ioVti.h"

#include <Teuchos_FileInputSource.hpp>
#include <EpetraExt_Utils.h> // for toString
#include <Teuchos_XMLParameterListWriter.hpp>


// =============================================================================
// Constructor
IoVti::IoVti( StaggeredGrid &sGrid ):
  SGrid(sGrid) // make this &sGrid
{
}
// =============================================================================



// =============================================================================
// Destructor
IoVti::~IoVti()
{
}
// =============================================================================



// =============================================================================
void IoVti::read( const std::string           &filename,
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
void IoVti::write( const std::string                 &filename,
                   const std::vector<double_complex> &psi,
                   const Teuchos::ParameterList      &problemParams )
{
  int    Nx = SGrid.getNx(),
         k,
         index[2];
  double h  = SGrid.getH();

  std::string str;

  Teuchos::XMLObject xmlPointData("PointData");
  xmlPointData.addAttribute( "Scalars", "abs(psi)" );

  // first build the XML structure
  Teuchos::XMLObject xmlDataArrayAbs("DataArray");
  xmlDataArrayAbs.addAttribute( "type", "Float32" );
  xmlDataArrayAbs.addAttribute( "Name", "abs(psi)" );
  xmlDataArrayAbs.addAttribute( "format", "ascii" );
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = SGrid.i2k( index );
          xmlDataArrayAbs.addContent( EpetraExt::toString(abs(psi[k])) + " ");
      }
  }
  xmlPointData.addChild(xmlDataArrayAbs);

  Teuchos::XMLObject xmlDataArrayArg("DataArray");
  xmlDataArrayArg.addAttribute( "type", "Float32" );
  xmlDataArrayArg.addAttribute( "Name", "arg(psi)" );
  xmlDataArrayArg.addAttribute( "format", "ascii" );
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = SGrid.i2k( index );
          xmlDataArrayArg.addContent( EpetraExt::toString(arg(psi[k])) + " " );
      }
  }
  xmlPointData.addChild(xmlDataArrayArg);

  Teuchos::XMLObject xmlPiece("Piece");
  str = "0 " + EpetraExt::toString(Nx) + " 0 " + EpetraExt::toString(Nx) + " 0 0";
  xmlPiece.addAttribute( "Extent", str );
  xmlPiece.addChild(xmlPointData);

  Teuchos::XMLObject xmlImageData("ImageData");
  xmlImageData.addAttribute( "WholeExtent", str );
  xmlImageData.addAttribute( "Origin", "0 0 0" );
  str = EpetraExt::toString(h) + " " + EpetraExt::toString(h) + " 0";
  xmlImageData.addAttribute( "Spacing", str );
  xmlImageData.addChild(xmlPiece);


  // append the problem parameters in XML form to the file
  Teuchos::XMLObject xmlParameterList("");
  xmlParameterList = Teuchos::XMLParameterListWriter::XMLParameterListWriter()
                                                        .toXML( problemParams );

  // define top level object
  Teuchos::XMLObject vtuxml("VTKFile");

  // append the parameter list to the embracing VTK XML object
  vtuxml.addChild(xmlParameterList);

  vtuxml.addAttribute( "type", "ImageData" );
  vtuxml.addAttribute( "version", "0.1" );
  vtuxml.addAttribute( "byte_order", "LittleEndian" );
  vtuxml.addChild(xmlImageData);

  // ---------------------------------------------------------------------------
  // write the contents to the file
  // open the file
  std::ofstream  vtkfile;
  vtkfile.open( filename.c_str() );

  // Do not plot the XML header as Teuchos' XML reader can't deal with it
  // vtkfile << "<?xml version=\"1.0\"?>" << std::endl;

  vtkfile << vtuxml;
  // close the file
  vtkfile.close();
  // ---------------------------------------------------------------------------
}
// =============================================================================




// =============================================================================
// Inside an XML object, this function looks for a specific tag and returns
// a pointer to it.
const Teuchos::XMLObject* IoVti::xmlFind ( const Teuchos::XMLObject *xmlObj,
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
const Teuchos::XMLObject* IoVti::xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
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