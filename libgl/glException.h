#ifndef GLEXCEPTION_H
#define GLEXCEPTION_H

// standard exceptions
#include <iostream>
#include <exception>
#include <string>

class glException: public std::exception
{

  public:

      //! Constructor.
      //! @param functionName Name of the function where the exception is thrown.
      //! @param errorMessage Error message for the exception.
      glException( const std::string& functionName,
                   const std::string& errorMessage  ) throw();

      //! Destructor.
      virtual ~glException() throw();

      //! Prints function name and error message to cout in a formatted way.
      virtual const char* what() const throw();

  protected:

  private:
      std::string _functionName;
      std::string _errorMessage;

};
#endif // GLEXCEPTION_H