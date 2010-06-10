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

#include "Ginla_StatusTest_ParameterLimits.h"

// #include "Ginla_Helpers.h"

#include <LOCA_Stepper.H>
#include <NOX_Utils.H>
// #include <NOX_Epetra_Group.H>
// #include <LOCA_MultiContinuation_AbstractGroup.H>


// ============================================================================
Ginla::StatusTest::ParameterLimits::
ParameterLimits( double lowerLimit,
                 double upperLimit,
                 bool return_failed_on_max_steps
               ) :
  lowerLimit_( lowerLimit ),
  upperLimit_( upperLimit ),
  tol_( 1.0e-10 ),
  value_( 0.0 ),
  return_failed_on_max_steps_( return_failed_on_max_steps ),
  status_( LOCA::StatusTest::Unevaluated )
{
}
// ============================================================================
Ginla::StatusTest::ParameterLimits::
~ParameterLimits()
{
}
// ============================================================================
LOCA::StatusTest::StatusType
Ginla::StatusTest::ParameterLimits::
checkStatus( const LOCA::Stepper& stepper,
                   LOCA::StatusTest::CheckType checkType )
{
  switch (checkType)
  {
  case LOCA::StatusTest::Complete:
  case LOCA::StatusTest::Minimal:
    value_ = stepper.getContinuationParameter();
    if ( stepper.getStepNumber() > 0 && // don't test on the first step
         (value_ < lowerLimit_ - tol_ || value_ > upperLimit_ + tol_) )
    {
      if ( return_failed_on_max_steps_ )
        status_ = LOCA::StatusTest::Failed;
      else
        status_ = LOCA::StatusTest::Finished;
    }
    else
    {
      status_ = LOCA::StatusTest::NotFinished;
    }
    break;
  case LOCA::StatusTest::None:
  default:
    status_ = LOCA::StatusTest::Unevaluated;
    break;
  }

  return status_;
}
// ============================================================================
LOCA::StatusTest::StatusType
Ginla::StatusTest::ParameterLimits::
getStatus() const
{
  return status_;
}
// ============================================================================
ostream&
Ginla::StatusTest::ParameterLimits::
print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status_;
  stream << NOX::Utils::sciformat(lowerLimit_,3)
         << " < "
         << "Parameter value = " << NOX::Utils::sciformat(value_,3)
         << " < "
         << NOX::Utils::sciformat(upperLimit_,3)
         << std::endl;
 return stream;
}
// ============================================================================