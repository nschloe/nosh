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

#include "Ginla_StatusTest_StabilityChange.h"

// #include "Ginla_StatusTest_Energy.h"
// 
// #include "Ginla_LocaSystem_Bordered.h"
// #include "Ginla_Helpers.h"
// 
// #include <LOCA_Stepper.H>
// #include <NOX_Epetra_Group.H>
// #include <LOCA_MultiContinuation_AbstractGroup.H>


// ============================================================================
Ginla::StatusTest::StabilityChange::
StabilityChange( const Teuchos::RCP<const Teuchos::ParameterList> & eigenInfo,
                 const int                                          stabilityChangeThreshold
               ):
    eigenInfo_( eigenInfo ),
    stabilityChangeThreshold_( stabilityChangeThreshold ),
    refNumUnstableEigenvalues_( 0 ),
    numUnstableEigenvalues_( 0 ),
    firstTime_( true ),
    status_( LOCA::StatusTest::Unevaluated )
{
}
// ============================================================================
Ginla::StatusTest::StabilityChange::
~StabilityChange()
{
}
// ============================================================================
LOCA::StatusTest::StatusType
Ginla::StatusTest::StabilityChange::
checkStatus( const LOCA::Stepper& stepper,
                   LOCA::StatusTest::CheckType checkType )
{
  // fetch the number of unstable eigenvalues
  numUnstableEigenvalues_ = eigenInfo_->get<unsigned int>( "#0unstable" );
  
  if ( firstTime_ )
  {
    refNumUnstableEigenvalues_ = numUnstableEigenvalues_;
    firstTime_ = false;
  }
  
  switch (checkType)
  {
  case LOCA::StatusTest::Complete:
  case LOCA::StatusTest::Minimal:

    if ( abs(numUnstableEigenvalues_-refNumUnstableEigenvalues_) >= stabilityChangeThreshold_ )
        status_ = LOCA::StatusTest::Finished;
    else
        status_ = LOCA::StatusTest::NotFinished;
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
Ginla::StatusTest::StabilityChange::
getStatus() const
{
  return status_;
}
// ============================================================================
ostream&
Ginla::StatusTest::StabilityChange::
print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
      stream << ' ';
  stream << status_;
  stream << "#unstable ev = " << numUnstableEigenvalues_;
  stream << std::endl;
  return stream;
}
// ============================================================================
