/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

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

#include <Teuchos_RCP.hpp>
#include <LOCA_StatusTest_Abstract.H>

// forward declarations
namespace Ginla {
  namespace LocaSystem {
    class Bordered;
  }
}
namespace Recti {
  namespace Grid {
    class General;
  }
}


namespace Ginla {

namespace StatusTest {

class Energy:
    public LOCA::StatusTest::Abstract
{
public:
  //! Constructor.
  Energy( const Teuchos::RCP<const Ginla::LocaSystem::Bordered> & glSystem,
          const Teuchos::RCP<const Recti::Grid::General>        & grid
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
  double tol_;
  double freeEnergy_;
  LOCA::StatusTest::StatusType status_;
  const Teuchos::RCP<const Ginla::LocaSystem::Bordered> glSystem_;
  const Teuchos::RCP<const Recti::Grid::General> grid_;

};

}

}

#endif // GINLA_STATUSTEST_ENERGY_H
