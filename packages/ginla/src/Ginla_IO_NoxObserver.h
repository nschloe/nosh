/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Sch\"omer

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

#ifndef GINLA_IO_NOXOBSERVER_H
#define GINLA_IO_NOXOBSERVER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Piro_Epetra_NOXObserver.hpp>

// forward declarations
namespace Ginla {
  class Komplex;
  namespace IO {
    class StateWriter;
  }
  namespace ModelEvaluator{
    class Default;
  }
}
namespace Recti {
  namespace Grid {
    class General;
  }
}


namespace Ginla {

namespace IO {

class NoxObserver:
    public Piro::Epetra::NOXObserver
{
public:
  //! Constructor
  NoxObserver ( const Teuchos::RCP<const Ginla::IO::StateWriter>         & stateWriter,
                const Teuchos::RCP<const Ginla::ModelEvaluator::Default> & modelEvaluator
              );
  
  //! Destructor
  virtual
  ~NoxObserver ();

  virtual
  void
  observeSolution(const Epetra_Vector& soln);
  
protected:
private:
    const Teuchos::RCP<const Ginla::IO::StateWriter>          stateWriter_;
    const Teuchos::RCP<const Ginla::ModelEvaluator::Default>  modelEvaluator_;
};

} // namespace IO
} // namespace Ginla

#endif // GINLA_IO_NOXOBSERVER_H