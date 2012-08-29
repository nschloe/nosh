// @HEADER
//
//    Helper class for writing statistics and states.
//    Copyright (C) 2010--2012  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
#ifndef NOSH_NOXOBSERVER_H
#define NOSH_NOXOBSERVER_H
// =============================================================================
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <NOX_Epetra_Observer.H>

#include "Nosh_CsvWriter.hpp"
// =============================================================================
// forward declarations
namespace Komplex2 {
class LinearProblem;
}
namespace Nosh {
namespace ModelEvaluator {
class Virtual;
}
}
// =============================================================================
namespace Nosh {

class Observer :
  public NOX::Epetra::Observer
{
public:
enum EObserverType { OBSERVER_TYPE_NEWTON,
                     OBSERVER_TYPE_CONTINUATION,
                     OBSERVER_TYPE_TURNING_POINT };

public:
//! Constructor
Observer(const Teuchos::RCP<const Nosh::ModelEvaluator::Virtual> &modelEval,
         const std::string & filename,
         const Observer::EObserverType &problemType
         );

//! Destructor
virtual
~Observer ();

virtual
void
observeSolution( const Epetra_Vector &soln );

protected:
private:
void
observeContinuation_(const Epetra_Vector &soln);

void
observeTurningPointContinuation_(const Epetra_Vector &soln);

void
saveContinuationStatistics_(const int stepIndex,
                            const Epetra_Vector &soln
                            );

private:

const Teuchos::RCP<const Nosh::ModelEvaluator::Virtual> modelEval_;
Nosh::CsvWriter csvWriter_;
const EObserverType observerType_;

};

} // namespace Nosh

#endif // NOSH_NOXOBSERVER_H
