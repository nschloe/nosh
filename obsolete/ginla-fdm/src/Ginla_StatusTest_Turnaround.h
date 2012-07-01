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

#ifndef GINLA_STATUSTEST_TURNAROUND_H
#define GINLA_STATUSTEST_TURNAROUND_H

#include "Ginla_Typedefs.h"

#include <Teuchos_RCP.hpp>
#include <LOCA_StatusTest_Abstract.H>

// forward declarations
namespace NOX {
  namespace Abstract {
    class Vector;
  }
}


namespace Ginla {

namespace StatusTest {

class Turnaround:
    public LOCA::StatusTest::Abstract
{
public:
  //! Constructor.
  Turnaround();

  //! Destructor.
  virtual
  ~Turnaround();

  virtual
  LOCA::StatusTest::StatusType
  checkStatus( const LOCA::Stepper& stepper,
                     LOCA::StatusTest::CheckType checkType
             );

  //! Return the result of the most recent checkStatus call
  virtual
  LOCA::StatusTest::StatusType
  getStatus() const;

  //! Output formatted description of stopping test to output stream.
  virtual
  ostream&
  print( ostream& stream,
         int indent = 0 ) const;

protected:

private:
  void
  computeUpdateProjection( const LOCA::Stepper & stepper );
  
private:
  bool firstTime_;
  double tol_;
  double updateProjection_;

  LOCA::StatusTest::StatusType status_;
  Teuchos::RCP<const NOX::Abstract::Vector> previousPoint_;
  Teuchos::RCP<NOX::Abstract::Vector> previousUpd_;

};

}

}

#endif // GINLA_STATUSTEST_TURNAROUND_H
