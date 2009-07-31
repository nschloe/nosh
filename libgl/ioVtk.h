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

  private:
      inline double_complex polar2complex( double abs,
                                           double arg  );

};