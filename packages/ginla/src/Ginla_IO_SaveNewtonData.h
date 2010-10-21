/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl\"omer

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

#ifndef GINLA_IO_SAVENEWTONDATA_H_
#define GINLA_IO_SAVENEWTONDATA_H_

#include <Teuchos_RCP.hpp>
#include <NOX_Abstract_PrePostOperator.H>

// foward declarations
namespace Ginla {
  namespace IO {
    class StateWriter;
  }
  namespace StateTranslator {
      class Virtual;
  }
}

namespace Ginla {
namespace IO {
      class SaveNewtonData:
        public NOX::Abstract::PrePostOperator
{

public:

  //! Constructor.
  SaveNewtonData ( const Teuchos::RCP<const Ginla::IO::StateWriter>          & stateWriter,
                   const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & translator
                 );

  //! Destructor.
  ~SaveNewtonData();

  //! Function that gets called before each iteration.
  //! This particular implementation prints the current state to the file
  //! data/newton-step-numRunPreIterate.vtk .
  //! @param solver The solver.
  void runPostIterate(const NOX::Solver::Generic& solver);

protected:
private:
  //! How ofter the function has been invoked yet.
  int numRunPreIterate;

  const Teuchos::RCP<const Ginla::IO::StateWriter>          stateWriter_;
  const Teuchos::RCP<const Ginla::StateTranslator::Virtual> translator_;
};

  } // namespace IO
} // namespace SaveNewtonData

#endif // GINLA_IO_SAVENEWTONDATA_H_
