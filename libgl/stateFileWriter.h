#ifndef STATEFILEWRITER_H
#define STATEFILEWRITER_H

#include "staggeredGrid.h"

#include <string>
#include <vector>
#include <complex>

#include <Teuchos_ParameterList.hpp>

typedef std::complex<double> double_complex;

class StateFileWriter
{
  public:

     //! Default constructor.
     StateFileWriter( StaggeredGrid &sGrid );

     //! Destructor
     ~StateFileWriter();

     //
     void writeFile( const std::string                 &fileName,
                     const std::string                 &fileFormat,
                     const std::vector<double_complex> &psi,
                     const Teuchos::ParameterList      &problemParams );

  private:

     StaggeredGrid::StaggeredGrid SGrid;

     void writeLegacyVtkFile( const std::string                 &filename,
                              const std::vector<double_complex> &psi,
                              const Teuchos::ParameterList      &problemParams );
     void writeVtiFile      ( const std::string                 &filename,
                              const std::vector<double_complex> &psi,
                              const Teuchos::ParameterList      &problemParams );
     void writeXdmfFile     ( const std::string                 &filename,
                              const std::vector<double_complex> &psi,
                              const Teuchos::ParameterList      &problemParams );

};
#endif // STATEFILEWRITER_H