#include "ioException.h"

// =============================================================================
IoException::IoException( const std::string & functionName,
                          const std::string & errorMessage  ) throw() :
  _functionName( functionName ),
  _errorMessage( errorMessage )
{
}
// =============================================================================
IoException::~IoException() throw()
{
}
// =============================================================================
const char* 
IoException::what() const throw()
{
  std::string message =  _functionName + ":\n"
                      + "\t" + _errorMessage + ".";

  return message.c_str();
}
// =============================================================================
