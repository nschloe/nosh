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

#include "Ginla_StatusTest_Energy.h"

#include "Ginla_State_Virtual.h"
#include "Ginla_StateTranslator.h"
#include "Ginla_Helpers.h"

#include <LOCA_Stepper.H>
#include <NOX_Epetra_Group.H>
#include <LOCA_MultiContinuation_AbstractGroup.H>


// ============================================================================
Ginla::StatusTest::Energy::
Energy( const Teuchos::RCP<const Ginla::StateTranslator> & stateTranslator,
        const double                                       maxFreeEnergy
      ) :
  freeEnergy_( 0.0 ),
  maxFreeEnergy_( maxFreeEnergy ),
  status_( LOCA::StatusTest::Unevaluated ),
  stateTranslator_( stateTranslator )
{
}
// ============================================================================
Ginla::StatusTest::Energy::
~Energy()
{
}
// ============================================================================
LOCA::StatusTest::StatusType
Ginla::StatusTest::Energy::
checkStatus( const LOCA::Stepper& stepper,
                   LOCA::StatusTest::CheckType checkType )
{
  
  switch (checkType)
  {
  case LOCA::StatusTest::Complete:
  case LOCA::StatusTest::Minimal:
    computeFreeEnergy( stepper );
    if ( freeEnergy_ > maxFreeEnergy_ )
      status_ = LOCA::StatusTest::Finished; // LOCA::StatusTest::Failed
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
Ginla::StatusTest::Energy::
computeFreeEnergy( const LOCA::Stepper & stepper )
{
    Teuchos::RCP<const NOX::Abstract::Group> solGroup =
        Teuchos::rcp_dynamic_cast<const NOX::Abstract::Group> ( stepper.getSolutionGroup() );
    const Epetra_Vector & x =
        ( Teuchos::dyn_cast<const NOX::Epetra::Vector> ( solGroup->getX() ) ).getEpetraVector();
        
    freeEnergy_ = stateTranslator_->createState( x )->freeEnergy();
    
    return;
}
// ============================================================================
LOCA::StatusTest::StatusType
Ginla::StatusTest::Energy::
getStatus() const
{
  return status_;
}
// ============================================================================
ostream&
Ginla::StatusTest::Energy::
print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
      stream << ' ';
  stream << status_;
  stream << "Free energy = " << NOX::Utils::sciformat( freeEnergy_, 3 );
  stream << " < " << NOX::Utils::sciformat(maxFreeEnergy_,3);
  stream << std::endl;
 return stream;
}
// ============================================================================