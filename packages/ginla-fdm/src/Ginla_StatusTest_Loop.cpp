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

#include "Ginla_StatusTest_Loop.h"

#include "Ginla_StatusTest_Energy.h"

// #include "Ginla_LocaSystem_Bordered.h"
#include "Ginla_Helpers.h"

#include <LOCA_Stepper.H>
#include <NOX_Epetra_Group.H>
#include <LOCA_MultiContinuation_AbstractGroup.H>


// ============================================================================
Ginla::StatusTest::Loop::
        Loop( const Teuchos::RCP<const Ginla::StateTranslator::Virtual> & stateTranslator ):
    firstTime_( true ),
    wasAway_( false ),
    tol_( 1.0e-3 ), // TODO make this depend on the maximum step size
    diffNorm_( 0.0 ),
    stateTranslator_( stateTranslator ),
    status_( LOCA::StatusTest::Unevaluated ),
    referenceState_( Teuchos::null )
{
}
// ============================================================================
Ginla::StatusTest::Loop::
~Loop()
{
}
// ============================================================================
LOCA::StatusTest::StatusType
Ginla::StatusTest::Loop::
checkStatus( const LOCA::Stepper               & stepper,
                   LOCA::StatusTest::CheckType   checkType )
{
  // store the reference solution
  if ( firstTime_ )
  {
      setReferencePoint( stepper );
      status_ = LOCA::StatusTest::NotFinished;
      firstTime_ = false;
      return status_;
  }

  switch (checkType)
  {
  case LOCA::StatusTest::Complete:
  case LOCA::StatusTest::Minimal:
    computeDiffNorm( stepper );
    if ( diffNorm_ > tol_ )
        wasAway_ = true;

    if ( !wasAway_ )
        status_ = LOCA::StatusTest::Unevaluated;
    else
    {
        if ( diffNorm_ < tol_ )
          status_ = LOCA::StatusTest::Finished;
        else
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
void
Ginla::StatusTest::Loop::
setReferencePoint( const LOCA::Stepper & stepper )
{
    Teuchos::RCP<const NOX::Abstract::Group> solGroup =
        Teuchos::rcp_dynamic_cast<const NOX::Abstract::Group> ( stepper.getSolutionGroup() );
    const Epetra_Vector & x =
        ( Teuchos::dyn_cast<const NOX::Epetra::Vector> ( solGroup->getX() ) ).getEpetraVector();

    referenceState_ = stateTranslator_->createState( x );

    return;
}
// ============================================================================
void
Ginla::StatusTest::Loop::
computeDiffNorm( const LOCA::Stepper & stepper )
{
    Teuchos::RCP<const NOX::Abstract::Group> solGroup =
        Teuchos::rcp_dynamic_cast<const NOX::Abstract::Group> ( stepper.getSolutionGroup() );
    const Epetra_Vector & x =
        ( Teuchos::dyn_cast<const NOX::Epetra::Vector> ( solGroup->getX() ) ).getEpetraVector();

    const Teuchos::RCP<Ginla::State::Updatable> state = stateTranslator_->createState( x );

    TEUCHOS_ASSERT( !referenceState_.is_null() );

    state->update( 1.0, *referenceState_, -1.0 );
    diffNorm_ = state->normalizedScaledL2Norm();

    return;
}
// ============================================================================
LOCA::StatusTest::StatusType
Ginla::StatusTest::Loop::
getStatus() const
{
  return status_;
}
// ============================================================================
ostream&
Ginla::StatusTest::Loop::
print( ostream& stream,
       int indent
     ) const
{
  for (int j = 0; j < indent; j ++)
      stream << ' ';

  stream << status_;
  stream << "||psi-ref||_L2 = " << NOX::Utils::sciformat(fabs(diffNorm_),3);
  stream << " > " << NOX::Utils::sciformat(tol_,3);
  stream << std::endl;
 return stream;
}
// ============================================================================
