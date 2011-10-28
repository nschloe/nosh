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

#ifndef GINLA_STATUSTEST_PARAMETERLIMITS_H
#define GINLA_STATUSTEST_PARAMETERLIMITS_H

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include "LOCA_StatusTest_Abstract.H"

#include "Teuchos_RCP.hpp"
#include "LOCA_GlobalData.H"

namespace Ginla {

namespace StatusTest {

class ParameterLimits: public LOCA::StatusTest::Abstract {

public:

  //! Constructor.
  ParameterLimits( double lowerLimit,
                   double upperLimit,
                   bool return_failed_on_max_steps = true
                  );

  //! Destructor.
  virtual ~ParameterLimits();

  virtual LOCA::StatusTest::StatusType
  checkStatus(const LOCA::Stepper& stepper,
              LOCA::StatusTest::CheckType status);

  virtual LOCA::StatusTest::StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  double lowerLimit_;
  double upperLimit_;

  double tol_;
  double value_;

  const bool return_failed_on_max_steps_;

  LOCA::StatusTest::StatusType status_;

  Teuchos::RCP<const LOCA::GlobalData> globalDataPtr_;

};

} // namespace Status
} // namespace LOCA

#endif // GINLA_STATUSTEST_PARAMETERLIMITS_H
