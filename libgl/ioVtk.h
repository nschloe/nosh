#include "ioVirtual.h"

class IoVtk: public IoVirtual
{
  public:

     //! Default constructor.
     IoVtk( std::string fname );

     //! Destructor
     virtual ~IoVtk();

     virtual void read( std::vector<double_complex> *psi,
                        Teuchos::ParameterList      *problemParams );

     virtual void write( const std::vector<double_complex> &psi,
                         const Teuchos::ParameterList      &problemParams,
                         StaggeredGrid                     &sGrid          );

  protected:
  private:

      //! joins a vector of strings to one string with a separator string sep
      std::string strJoin( const std::vector<std::string> & vec,
                           const std::string              & sep  );

};