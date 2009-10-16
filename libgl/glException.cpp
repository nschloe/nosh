#include "glException.h"

// =============================================================================
glException::glException( const std::string & functionName,
                          const std::string & errorMessage  ) throw() :
  _functionName( functionName ),
  _errorMessage( errorMessage )
{
}
// =============================================================================
glException::~glException() throw()
{
}
// =============================================================================
const char* glException::what() const throw()
{
  std::string message =  _functionName + ":\n"
                      + "\t" + _errorMessage + ".";

  return message.c_str();
}
// =============================================================================