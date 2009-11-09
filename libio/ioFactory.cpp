#include "ioFactory.h"
#include "ioException.h"

#include "ioVtk.h"
#include "ioVti.h"
#include "ioXdmf.h"

// =============================================================================
IoVirtual*
IoFactory::createFileIo( std::string fileName )
{

  // analyze the file name for extension
  int         dotPos    = fileName.rfind(".");
  std::string extension = fileName.substr( dotPos+1, fileName.size()-dotPos-1  );

  if ( extension.compare("vtk") == 0) {
      return new IoVtk( fileName );
  } else if (extension.compare("vti") == 0) {
      return new IoVti( fileName );
  } else if (extension.compare("xmf") == 0) {
      return new IoXdmf( fileName );
  } else {
      std::string message = "File name extension \"" + extension + "\" "
                          + "not recognized. Must be one of \"vtk\", "
                          + "\"vti\", \"xmf\".";
      throw IoException( "IoFactory::createFileIo",
                          message );
  }

}
// =============================================================================
