/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef GINLA_BRANCHEXPLORER_DRIVER_H
#define GINLA_BRANCHEXPLORER_DRIVER_H

#include <string>

#include "Ginla_Typedefs.h"

#include <Epetra_Comm.h>

// forward declarations
namespace Recti {
  namespace Grid {
    class Uniform;
  }
}

namespace Ginla {
namespace BranchExplorer {

class Driver
{
public:
  //! Constructor.
  Driver();
  
  //! Destructor.
  virtual
  ~Driver();
  
  //! Start the branch explorer.
  void
  run();
  
protected:
  
private:
  void
  oneContinuation( const Teuchos::RCP<Epetra_Comm>            & eComm,
                   const std::string                          & outputDirectory,
                   const Teuchos::RCP<const ComplexVector>    & psi,
                   const Teuchos::RCP<Recti::Grid::Uniform>   & grid,
                   const Teuchos::ParameterList               & glParameters,
                   const Teuchos::RCP<Teuchos::ParameterList> & paramList
                 );
  
private:
   const std::string contFileBaseName_;
   const std::string outputFormat_;
   const std::string contDataFileName_;
};

}
}

#endif // GINLA_BRANCHEXPLORER_DRIVER_H
