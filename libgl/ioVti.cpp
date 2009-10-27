#include "ioVti.h"

#include "glException.h"

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
void IoVti::read( Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > &psi,
                  const Teuchos::RCP<const Teuchos::Comm<int> >          comm, // TODO: remove this
                  Teuchos::ParameterList                                 &problemParams
                ) const
{

  throw glException( "IoVti::StateFileReader",
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

    if ( !absPsiObject ) { // pointer invalid
        throw glException( "GlSystem::GlSystem",
                           "No such XML Object found." );
    }

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
void IoVti::write( const Tpetra::MultiVector<double_complex,int> &psi,
                   const Teuchos::ParameterList                  &problemParams,
                   const StaggeredGrid                           &sGrid
                 ) const
{
  int    Nx = sGrid.getNx();
  int    k;
  double h  = sGrid.getH();

  std::string str;

  Teuchos::XMLObject xmlPointData("PointData");
  xmlPointData.addAttribute( "Scalars", "abs(psi)" );

  // first build the XML structure
  Teuchos::XMLObject xmlDataArrayAbs("DataArray");
  xmlDataArrayAbs.addAttribute( "type", "Float32" );
  xmlDataArrayAbs.addAttribute( "Name", "abs(psi)" );
  xmlDataArrayAbs.addAttribute( "format", "ascii" );
  
  Teuchos::ArrayRCP<const double_complex> psiView = psi.getVector(0)->get1dView();
  Teuchos::Array<int> index(2);
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = sGrid.i2k( index );
          xmlDataArrayAbs.addContent( EpetraExt::toString(abs(psiView[k])) + " ");
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
          k = sGrid.i2k( index );
          xmlDataArrayArg.addContent( EpetraExt::toString(arg(psiView[k])) + " " );
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
IoVti::write( const Tpetra::MultiVector<double_complex,int> &psi,
              const StaggeredGrid                           &sGrid
            ) const
{
    throw glException( "IoVti::write", "Method not yet implemented." );
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
