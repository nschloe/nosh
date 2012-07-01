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
#ifndef CUANTICO_NOXOBSERVER_H
#define CUANTICO_NOXOBSERVER_H

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <NOX_Epetra_Observer.H>
// =============================================================================
// forward declarations
namespace Komplex2 {
class LinearProblem;
}
namespace Cuantico {
class CsvWriter;
class ModelEvaluator;
class State;
}
// =============================================================================
namespace Cuantico {

class NoxObserver :
  public NOX::Epetra::Observer
{
public:
enum EObserverType { OBSERVER_TYPE_NEWTON,
                     OBSERVER_TYPE_CONTINUATION,
                     OBSERVER_TYPE_TURNING_POINT };

public:
//! Constructor
NoxObserver(const Teuchos::RCP<const Cuantico::ModelEvaluator> &modelEval,
            const std::string & filename,
            const NoxObserver::EObserverType &problemType
            );

//! Destructor
virtual
~NoxObserver ();

virtual
void
observeSolution( const Epetra_Vector &soln );

protected:
private:
void
observeContinuation_( const Teuchos::RCP<const Cuantico::State> &state
                      );
void
observeTurningPointContinuation_( const Teuchos::RCP<const Cuantico::State> &state
                                  );

void
saveContinuationStatistics_( const int stepIndex,
                             const Teuchos::RCP<const Cuantico::State> &state
                             );

private:
private:

const Teuchos::RCP<const Cuantico::ModelEvaluator> modelEval_;
const Teuchos::RCP<Cuantico::CsvWriter> csvWriter_;
const EObserverType observerType_;

};

} // namespace Cuantico

#endif // CUANTICO_NOXOBSERVER_H
