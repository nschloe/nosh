#include "staggeredGrid.h"

#include <string>
#include <vector>
#include <complex>

#include <Teuchos_ParameterList.hpp>

typedef std::complex<double> double_complex;

class IoXdmf
{
  public:

     //! Default constructor.
     IoXdmf( StaggeredGrid &sGrid );

     //! Destructor
     ~IoXdmf();

     void read( const std::string           &fileName,
                std::vector<double_complex> *psi,
                Teuchos::ParameterList      *problemParams );

     //
     void write( const std::string                 &fileName,
                 const std::vector<double_complex> &psi,
                 const Teuchos::ParameterList      &problemParams );

  private:

     StaggeredGrid::StaggeredGrid SGrid;

};