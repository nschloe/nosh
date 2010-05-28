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

#include "Ginla_StatusTest_MaxAcceptedSteps.h"

// #include "LOCA_StatusTest_Abstract.H"

#include "LOCA_Stepper.H"

// FIXME remove these headers?
#include "NOX_Utils.H"
#include "LOCA_GlobalData.H"

Ginla::StatusTest::MaxAcceptedSteps::
MaxAcceptedSteps( int maxAcceptedSteps,
                  bool return_failed_on_max_steps,
                  const Teuchos::RCP<const LOCA::GlobalData> globalDataPtr
                ) :
  maxAcceptedSteps_(maxAcceptedSteps),
  return_failed_on_max_steps_(return_failed_on_max_steps),
  numAcceptedSteps_(0),
  status(LOCA::StatusTest::Unevaluated)
{
  if ( globalDataPtr.is_valid_ptr() && !globalDataPtr.is_null() )
    globalDataPtr_ = globalDataPtr;

  if (maxAcceptedSteps_ < 0)
  {
    if ( globalDataPtr_.is_valid_ptr() && !globalDataPtr_.is_null() )
        globalDataPtr_->locaUtils->err() << "LOCA::StatusTest::MaxIters - must choose a number greater than or equal to zero" << endl;
    else
        // This will spit out the error message NUMPROC times. -- Without locaUtils, there's nothing we can do..
        std::cerr << "LOCA::StatusTest::MaxIters - must choose a number greater than or equal to zero" << endl;
    throw "LOCA Error";
  }
}

Ginla::StatusTest::MaxAcceptedSteps::
~MaxAcceptedSteps()
{
}

LOCA::StatusTest::StatusType Ginla::StatusTest::MaxAcceptedSteps::
checkStatus(const LOCA::Stepper& stepper,
            LOCA::StatusTest::CheckType checkType)
{
  switch (checkType)
  {
  case LOCA::StatusTest::Complete:
  case LOCA::StatusTest::Minimal:
    numAcceptedSteps_ = stepper.getStepNumber();
    if (numAcceptedSteps_ >= maxAcceptedSteps_)
      status = return_failed_on_max_steps_ ? LOCA::StatusTest::Failed : LOCA::StatusTest::Finished;
    else
      status = LOCA::StatusTest::NotFinished;
    break;

  case LOCA::StatusTest::None:
  default:
    numAcceptedSteps_ = -1;
    status = LOCA::StatusTest::Unevaluated;
    break;
  }

  return status;
}

LOCA::StatusTest::StatusType Ginla::StatusTest::MaxAcceptedSteps::
getStatus() const
{
  return status;
}

ostream& Ginla::StatusTest::MaxAcceptedSteps::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Number of Iterations = " << numAcceptedSteps_ << " < " << maxAcceptedSteps_;
  stream << endl;
 return stream;
}

int Ginla::StatusTest::MaxAcceptedSteps::
getMaxAcceptedSteps() const
{
  return maxAcceptedSteps_;
}

int Ginla::StatusTest::MaxAcceptedSteps::
getNumAcceptedSteps() const
{
  return numAcceptedSteps_;
}
