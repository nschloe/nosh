#ifndef IOEXCEPTION_H
#define IOEXCEPTION_H

// standard exceptions
#include <iostream>
#include <exception>
#include <string>

class IoException: public std::exception
{

  public:

      //! Constructor.
      //! @param functionName Name of the function where the exception is thrown.
      //! @param errorMessage Error message for the exception.
	  IoException( const std::string& functionName,
                   const std::string& errorMessage  ) throw();

      //! Destructor.
      virtual ~IoException() throw();

      //! Prints function name and error message to cout in a formatted way.
      virtual const char* what() const throw();

  protected:

  private:
      std::string _functionName;
      std::string _errorMessage;

};
#endif // IOEXCEPTION_H