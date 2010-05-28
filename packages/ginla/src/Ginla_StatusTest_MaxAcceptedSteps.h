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

#ifndef GINLA_STATUSTEST_MAXACCEPTEDSTEPS_H
#define GINLA_STATUSTEST_MAXACCEPTEDSTEPS_H

#include "LOCA_StatusTest_Abstract.H"

#include "Teuchos_RCP.hpp"
#include "LOCA_GlobalData.H"

namespace Ginla {

namespace StatusTest {

class MaxAcceptedSteps: public LOCA::StatusTest::Abstract {

public:

  //! Constructor.
  MaxAcceptedSteps( int maxAcceptedSteps,
                    bool return_failed_on_max_steps = true,
                    const Teuchos::RCP<const LOCA::GlobalData> globalData = Teuchos::null
                  );

  //! Destructor.
  virtual ~MaxAcceptedSteps();

  virtual LOCA::StatusTest::StatusType
  checkStatus(const LOCA::Stepper& stepper,
              LOCA::StatusTest::CheckType status);

  virtual LOCA::StatusTest::StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

  virtual int getMaxAcceptedSteps() const;

  virtual int getNumAcceptedSteps() const;

private:

  int maxAcceptedSteps_;

  bool return_failed_on_max_steps_;

  int numAcceptedSteps_;

  LOCA::StatusTest::StatusType status;

  Teuchos::RCP<const LOCA::GlobalData> globalDataPtr_;

};

} // namespace Status
} // namespace LOCA

#endif // GINLA_STATUSTEST_MAXACCEPTEDSTEPS_H
