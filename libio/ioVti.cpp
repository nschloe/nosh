#include "ioVti.h"

#include <Teuchos_FileInputSource.hpp>
#include <EpetraExt_Utils.h> // for toString
#include <Teuchos_XMLParameterListWriter.hpp>


// =============================================================================
// Constructor
IoVti::IoVti( std::string fname ):
  IoVirtual(fname)
{
}
// =============================================================================
// Destructor
IoVti::~IoVti()
{
}
// =============================================================================
void
IoVti::read( const Teuchos::RCP<const Teuchos::Comm<int> >        &tComm,
	               Teuchos::RCP<Tpetra::MultiVector<double,int> > &x,
			       Teuchos::ParameterList                         &problemParams
	       ) const
{

  TEST_FOR_EXCEPTION( true,
	                  std::logic_error,
	                  "readVtiFile not yet implemented." );

  // pass a possible 
  //<?xml version="1.0"?>
  // at the beginning of the file


    Teuchos::FileInputSource xmlFile(fileName_);

    // Extract the object from the filename.
    // -- This is actually quite costly, as it read all -- *all* -- data in as
    // strings, including the heavy numerical data.
    Teuchos::XMLObject xmlFileObject = xmlFile.getObject();

//     // find and read the parameter list
//     const Teuchos::XMLObject* parameterListObject =
//                                       xmlFind( &xmlFileObject, "ParameterList" );
//     problemParams = Teuchos::XMLParameterListReader().
//                                           toParameterList( *parameterListObject );

    // find and get the absolute value vector
    std::string absPsiString;
    const Teuchos::XMLObject* absPsiObject = xmlAttributeFind( &xmlFileObject,
                                                              "DataArray",
                                                              "Name",
                                                              "abs(psi)"      );

    TEST_FOR_EXCEPTION( !absPsiObject, // pointer invalid
                        std::logic_error,
                        "No such XML Object found." );

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

}
// =============================================================================
void
IoVti::write ( const Tpetra::MultiVector<double,int> & x,
               const int                               Nx,
               const double                            h,
               const Teuchos::ParameterList          & problemParams
             ) const
{
  std::string str;

  Teuchos::XMLObject xmlPointData("PointData");
  xmlPointData.addAttribute( "Scalars", "abs(psi)" );

  int numVectors = x.getNumVectors();
  for (int k; k<numVectors; k++) {
	  // first build the XML structure
	  Teuchos::XMLObject xmlDataArray("DataArray");
	  xmlDataArray.addAttribute( "type", "Float32" );
	  std::string name = "x" + EpetraExt::toString(k);
	  xmlDataArray.addAttribute( "Name", name );
	  xmlDataArray.addAttribute( "format", "ascii" );

	  Teuchos::ArrayRCP<const double> xKView = x.getVector(k)->get1dView();
	  int l = 0;
	  for (int i=0; i<Nx+1; i++) {
		  for (int j=0; j<Nx+1; j++) {
			  l++;
			  xmlDataArray.addContent( EpetraExt::toString(abs(xKView[l])) + " ");
		  }
	  }
	  xmlPointData.addChild(xmlDataArray);
  }

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
  vtkfile.open( fileName_.c_str() );

  // Do not plot the XML header as Teuchos' XML reader can't deal with it
  // vtkfile << "<?xml version=\"1.0\"?>" << std::endl;

  vtkfile << vtuxml;
  // close the file
  vtkfile.close();
  // ---------------------------------------------------------------------------
}
// =============================================================================
void
IoVti::write( const Tpetra::MultiVector<double,int> & x,
              const int                               Nx,
              const double                            h
            ) const
{
  TEST_FOR_EXCEPTION( true,
	                  std::logic_error,
	                  "Not yet implemented." );
}
// =============================================================================
// Inside an XML object, this function looks for a specific tag and returns
// a pointer to it.
const Teuchos::XMLObject*
IoVti::xmlFind ( const Teuchos::XMLObject *xmlObj,
                 const std::string        tag
               ) const
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
const Teuchos::XMLObject*
IoVti::xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
                          const std::string        tag,
                          const std::string        attribute,
                          const std::string        value
                        ) const
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
