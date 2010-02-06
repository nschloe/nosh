#include "ioVti.h"

#include <Teuchos_FileInputSource.hpp>
#include <EpetraExt_Utils.h> // for toString
#include <Teuchos_XMLParameterListWriter.hpp>

#include <boost/filesystem/fstream.hpp>

// =============================================================================
// Constructor
IoVti::IoVti ( boost::filesystem::path fname ) :
        IoVirtual ( fname )
{
}
// =============================================================================
// Destructor
IoVti::~IoVti()
{
}
// =============================================================================
void
IoVti::read ( const Teuchos::RCP<const Teuchos::Comm<int> > &tComm,
              Teuchos::RCP<DoubleMultiVector> &x,
              Teuchos::ParameterList &problemParams ) const
{

    // TODO implement
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "readVtiFile not yet implemented." );

    // pass a possible
    //<?xml version="1.0"?>
    // at the beginning of the file

    Teuchos::FileInputSource xmlFile ( fileName_.string() );

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
    const Teuchos::XMLObject* absPsiObject = xmlAttributeFind ( &xmlFileObject,
            "DataArray", "Name", "abs(psi)" );

    TEST_FOR_EXCEPTION ( !absPsiObject, // pointer invalid
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
IoVti::read ( const Teuchos::RCP<const Teuchos::Comm<int> > &tComm,
              Teuchos::RCP<ComplexMultiVector> &x,
              Teuchos::ParameterList &problemParams ) const
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
void
IoVti::write ( const Epetra_MultiVector              & x,
               const Teuchos::Tuple<unsigned int,2>  & Nx,
               const Teuchos::Tuple<double,2>        & h,
               const Teuchos::ParameterList          & problemParams
             )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
void
IoVti::write ( const DoubleMultiVector               & x,
               const Teuchos::Tuple<unsigned int,2>  & Nx,
               const Teuchos::Tuple<double,2>        & h,
               const Teuchos::ParameterList          & problemParams
             )
{
    std::string str;

    Teuchos::XMLObject xmlPointData ( "PointData" );
    xmlPointData.addAttribute ( "Scalars", "abs(psi)" );

    int numVectors = x.getNumVectors();
    for ( int k; k < numVectors; k++ )
    {
        // first build the XML structure
        Teuchos::XMLObject xmlDataArray ( "DataArray" );
        xmlDataArray.addAttribute ( "type", "Float32" );
        std::string name = "x" + EpetraExt::toString ( k );
        xmlDataArray.addAttribute ( "Name", name );
        xmlDataArray.addAttribute ( "format", "ascii" );

        Teuchos::ArrayRCP<const double> xKView = x.getVector ( k )->get1dView();
        int l = 0;
        for ( unsigned int i = 0; i < Nx[0] + 1; i++ )
        {
            for ( unsigned int j = 0; j < Nx[1] + 1; j++ )
            {
                l++;
                xmlDataArray.addContent ( EpetraExt::toString ( abs ( xKView[l] ) ) + " " );
            }
        }
        xmlPointData.addChild ( xmlDataArray );
    }

    Teuchos::XMLObject xmlPiece ( "Piece" );
    str = "0 " + EpetraExt::toString ( Nx[0] ) + " 0 " + EpetraExt::toString ( Nx[1] )
          + " 0 0";
    xmlPiece.addAttribute ( "Extent", str );
    xmlPiece.addChild ( xmlPointData );

    Teuchos::XMLObject xmlImageData ( "ImageData" );
    xmlImageData.addAttribute ( "WholeExtent", str );
    xmlImageData.addAttribute ( "Origin", "0 0 0" );
    str = EpetraExt::toString ( h[0] ) + " " + EpetraExt::toString ( h[1] ) + " 0";
    xmlImageData.addAttribute ( "Spacing", str );
    xmlImageData.addChild ( xmlPiece );

    // append the problem parameters in XML form to the file
    Teuchos::XMLObject xmlParameterList ( "" );
    xmlParameterList
    = Teuchos::XMLParameterListWriter::XMLParameterListWriter() .toXML (
          problemParams );

    // define top level object
    Teuchos::XMLObject vtuxml ( "VTKFile" );

    // append the parameter list to the embracing VTK XML object
    vtuxml.addChild ( xmlParameterList );

    vtuxml.addAttribute ( "type", "ImageData" );
    vtuxml.addAttribute ( "version", "0.1" );
    vtuxml.addAttribute ( "byte_order", "LittleEndian" );
    vtuxml.addChild ( xmlImageData );

    // ---------------------------------------------------------------------------
    // write the contents to the file
    // open the file
    boost::filesystem::ofstream vtkfile ( fileName_ );

    // Do not plot the XML header as Teuchos' XML reader can't deal with it
    // vtkfile << "<?xml version=\"1.0\"?>" << std::endl;

    vtkfile << vtuxml;
    // close the file
    vtkfile.close();
    // ---------------------------------------------------------------------------
}
// =============================================================================
void
IoVti::write ( const ComplexMultiVector               & x,
               const Teuchos::Tuple<unsigned int,2>  & Nx,
               const Teuchos::Tuple<double,2>        & h,
               const Teuchos::ParameterList          & problemParams
             )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
void
IoVti::write ( const DoubleMultiVector              & x,
               const Teuchos::Tuple<unsigned int,2> & Nx,
               const Teuchos::Tuple<double,2>       & h
             )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
void
IoVti::write ( const ComplexMultiVector              & x,
               const Teuchos::Tuple<unsigned int,2>  & Nx,
               const Teuchos::Tuple<double,2>        & h
             )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
// Inside an XML object, this function looks for a specific tag and returns
// a pointer to it.
const Teuchos::XMLObject*
IoVti::xmlFind ( const Teuchos::XMLObject *xmlObj, const std::string tag ) const
{
    const Teuchos::XMLObject* xmlOut = NULL;

    if ( !xmlObj->getTag().compare ( tag ) ) // strings are equal
        xmlOut = xmlObj;
    else
        for ( int k = 0; k < xmlObj->numChildren(); k++ )
        {
            xmlOut = xmlFind ( & ( xmlObj->getChild ( k ) ), tag ); // recursive call
            if ( xmlOut )
                break; // not the null pointer => return
        }

    return xmlOut;
}
// =============================================================================
const Teuchos::XMLObject*
IoVti::xmlAttributeFind ( const Teuchos::XMLObject *xmlObj,
                          const std::string tag, const std::string attribute, const std::string value ) const
{
    const Teuchos::XMLObject* xmlOut = NULL;

    // std::cout << "Looking for: " << tag << " " << attribute << " " <<  value << std::endl;
    // std::cout << "Found:       " << xmlObj->getTag() << " " << xmlObj->hasAttribute(attribute) << std::endl;
    // if ( xmlObj->hasAttribute(attribute) )
    //      std::cout << "1" << std::endl;
    //      std::cout << xmlObj->getAttribute(attribute) << std::endl;

    if ( !xmlObj->getTag().compare ( tag ) && xmlObj->hasAttribute ( attribute )
         && xmlObj->getAttribute ( attribute ).compare ( value ) )
        xmlOut = xmlObj;
    else
        for ( int k = 0; k < xmlObj->numChildren(); k++ )
        {
            xmlOut = xmlAttributeFind ( & ( xmlObj->getChild ( k ) ), // recursive call
                                        tag, attribute, value );
            if ( xmlOut )
                break; // not the null pointer => return
        }

    return xmlOut;
}
// =============================================================================
void
IoVti::write ( const Epetra_MultiVector                        & x,
               const Teuchos::Array<Teuchos::Tuple<double,2> > & loc,
               const Teuchos::ParameterList                    & problemParams
             )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
void
IoVti::write ( const Epetra_MultiVector             & x,
               const Teuchos::Tuple<unsigned int,2> & Nx,
               const Teuchos::Tuple<double,2>       & h,
               const Teuchos::Array<int>            & kBoundingBox,
               const Teuchos::ParameterList         & problemParams,
               const double                         & dummyValue
             )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Not yet implemented." );
}
// =============================================================================
