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

#ifndef GINLA_STATUSTEST_STABILITYCHANGE_H
#define GINLA_STATUSTEST_STABILITYCHANGE_H

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

class StabilityChange:
    public LOCA::StatusTest::Abstract
{
public:
  //! Constructor.
  StabilityChange( const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo,
                   const int                                          stabilityChangeThreshold
                 );

  //! Destructor.
  virtual
  ~StabilityChange();

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
         int indent = 0 ) const;

protected:

private:
  
private:
  const Teuchos::RCP<const Teuchos::ParameterList> eigenInfo_;
  const int stabilityChangeThreshold_;
  int refNumUnstableEigenvalues_;
  int numUnstableEigenvalues_;
  bool firstTime_;
  LOCA::StatusTest::StatusType status_;

};

}

}

#endif // GINLA_STATUSTEST_STABILITYCHANGE_H
