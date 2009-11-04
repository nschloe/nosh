#include "ioException.h"

// =============================================================================
ioException::ioException( const std::string & functionName,
                          const std::string & errorMessage  ) throw() :
  _functionName( functionName ),
  _errorMessage( errorMessage )
{
}
// =============================================================================
ioException::~ioException() throw()
{
}
// =============================================================================
const char* 
ioException::what() const throw()
{
  std::string message =  _functionName + ":\n"
                      + "\t" + _errorMessage + ".";

  return message.c_str();
}
// =============================================================================
