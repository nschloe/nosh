// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2010, 2011  Nico Schl\"omer
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
#ifndef GINLA_NOXOBSERVER_H
#define GINLA_NOXOBSERVER_H

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
namespace Ginla {
  class StateWriter;
  class StatsWriter;
  class ModelEvaluator;
  class State;
}
// =============================================================================
namespace Ginla {

class NoxObserver:
    public NOX::Epetra::Observer
{
public:
    enum EObserverType { OBSERVER_TYPE_NEWTON,
                         OBSERVER_TYPE_CONTINUATION,
                         OBSERVER_TYPE_TURNING_POINT };

public:
  //! Constructor
  NoxObserver ( const Teuchos::RCP<const Ginla::ModelEvaluator> & modelEval,
                const NoxObserver::EObserverType                & problemType
              );

  //! Destructor
  virtual
  ~NoxObserver ();

  virtual
  void
  observeSolution(const Epetra_Vector& soln);

  void
  setStatisticsWriter( const Teuchos::RCP<Ginla::StatsWriter> & statsWriter );

protected:
private:
  void
  observeContinuation_( const Teuchos::RCP<const Ginla::State> & state
                      );
  void
  observeTurningPointContinuation_( const Teuchos::RCP<const Ginla::State> & state
                                  );

  void
  saveContinuationStatistics_( const int stepIndex,
                               const Teuchos::RCP<const Ginla::State> & state
                             );

private:
private:

    const Teuchos::RCP<const Ginla::ModelEvaluator> modelEval_;
    const EObserverType observerType_;

    Teuchos::RCP<Ginla::StatsWriter> statsWriter_;
};

} // namespace Ginla

#endif // GINLA_NOXOBSERVER_H
