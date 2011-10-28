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

#ifndef GINLA_STATUSTEST_ENERGY_H
#define GINLA_STATUSTEST_ENERGY_H

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_RCP.hpp>
#include <LOCA_StatusTest_Abstract.H>

// forward declarations
namespace Ginla {
    namespace StateTranslator {
        class Virtual;
    }
}


namespace Ginla {

namespace StatusTest {

class Energy:
    public LOCA::StatusTest::Abstract
{
public:
  //! Constructor.
    Energy( const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & stateTranslator,
            const double                                                maxFreeEnergy
        );

  //! Destructor.
  virtual
  ~Energy();

  virtual
  LOCA::StatusTest::StatusType
  checkStatus(const LOCA::Stepper& stepper,
                    LOCA::StatusTest::CheckType checkType);

  //! Return the result of the most recent checkStatus call
  virtual
  LOCA::StatusTest::StatusType
  getStatus() const;

  //! Output formatted description of stopping test to output stream.
  virtual
  ostream&
  print( ostream& stream,
         int indent = 0) const;

protected:

private:
  void
  computeFreeEnergy( const LOCA::Stepper & stepper );

private:
  double freeEnergy_;
  double maxFreeEnergy_;
  LOCA::StatusTest::StatusType status_;
  const Teuchos::RCP<const Ginla::StateTranslator::Virtual> stateTranslator_;

};

}

}

#endif // GINLA_STATUSTEST_ENERGY_H
