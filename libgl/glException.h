// standard exceptions
#include <iostream>
#include <exception>
#include <string>

class glException: public std::exception
{

  public:

      glException( const std::string& functionName,
                   const std::string& errorMessage  ) throw();

      virtual ~glException() throw();

      virtual const char* what() const throw();

  protected:

  private:
      std::string _functionName;
      std::string _errorMessage;

};
