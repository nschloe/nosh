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

#include "Ginla_StatusTest_Turnaround.h"

#include "Ginla_StatusTest_Energy.h"

// #include "Ginla_LocaSystem_Bordered.h"
#include "Ginla_Helpers.h"

#include <LOCA_Stepper.H>
#include <NOX_Epetra_Group.H>
#include <LOCA_MultiContinuation_AbstractGroup.H>


// ============================================================================
Ginla::StatusTest::Turnaround::
Turnaround():
    firstTime_ ( true ),
    tol_( 1.0e-12 ),
    updateProjection_( 0.0 ),
    status_( LOCA::StatusTest::Unevaluated ),
    previousPoint_( Teuchos::null ),
    previousUpd_( Teuchos::null )
{
}
// ============================================================================
Ginla::StatusTest::Turnaround::
~Turnaround()
{
}
// ============================================================================
LOCA::StatusTest::StatusType
Ginla::StatusTest::Turnaround::
checkStatus( const LOCA::Stepper& stepper,
                   LOCA::StatusTest::CheckType checkType )
{
  // store the reference solution
  if ( firstTime_ )
  {
      Teuchos::RCP<const NOX::Abstract::Group> solGroup =
        Teuchos::rcp_dynamic_cast<const NOX::Abstract::Group> ( stepper.getSolutionGroup() );
      const NOX::Abstract::Vector & x = solGroup->getX();
      previousPoint_ = x.clone();
      
      // set the initial update to 0
      previousUpd_ = x.clone();
      previousUpd_->scale( 0.0 );
      
      firstTime_ = false;
      return LOCA::StatusTest::Unevaluated;
  }
  
  switch (checkType)
  {
  case LOCA::StatusTest::Complete:
  case LOCA::StatusTest::Minimal:
    computeUpdateProjection( stepper );
    if ( abs(updateProjection_+1.0) < tol_ )
      status_ = LOCA::StatusTest::Failed; // LOCA::StatusTest::Failed
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
void
Ginla::StatusTest::Turnaround::
computeUpdateProjection( const LOCA::Stepper & stepper )
{
    Teuchos::RCP<const NOX::Abstract::Group> solGroup =
        Teuchos::rcp_dynamic_cast<const NOX::Abstract::Group> ( stepper.getSolutionGroup() );
    const NOX::Abstract::Vector & x = solGroup->getX();
        
    // calculate update
    Teuchos::RCP<NOX::Abstract::Vector> upd =  x.clone();
    upd->update( 1.0, *previousPoint_, -1.0 );

    // normalize
    double alpha = upd->norm();
    if ( fabs(alpha) > 1.0e-15 )
        upd->scale( 1.0/alpha );

    // project update to previous update
    updateProjection_ = upd->innerProduct( *previousUpd_ );
//     upd->update( 1.0, *previousUpd_, -1.0 );
//     updateProjection_ = upd->norm();

    // save for next step
    previousUpd_ = upd;
    previousPoint_ = x.clone();
    
    return;
}
// ============================================================================
LOCA::StatusTest::StatusType
Ginla::StatusTest::Turnaround::
getStatus() const
{
  return status_;
}
// ============================================================================
ostream&
Ginla::StatusTest::Turnaround::
print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
      stream << ' ';

  stream << status_;
  stream << "<upd,prevUpd>_2 = " << NOX::Utils::sciformat(updateProjection_,3);
  stream << std::endl;
 return stream;
}
// ============================================================================
