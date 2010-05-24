/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Sch\"omer

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

#ifndef GINLA_IO_NOXOBSERVER_H
#define GINLA_IO_NOXOBSERVER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <NOX_Epetra_Observer.H>

// forward declarations
namespace Ginla {
  class State;
  namespace Komplex {
    class LinearProblem;
  }
  namespace IO {
    class StateWriter;
    class StatsWriter;
  }
  namespace ModelEvaluator{
    class Default;
  }
  namespace Operator {
    class Virtual;
  }
}
namespace Recti {
  namespace Grid {
    class General;
  }
}


namespace Ginla {

namespace IO {

class NoxObserver:
    public NOX::Epetra::Observer
{
public:
    enum ObserverType { NONLINEAR,
                        CONTINUATION,
                        TURNING_POINT };

public:
  //! Constructor
  NoxObserver ( const Teuchos::RCP<const Ginla::IO::StateWriter>         & stateWriter,
                const Teuchos::RCP<const Ginla::ModelEvaluator::Default> & modelEvaluator,
                const NoxObserver::ObserverType                          & problemType
              );
  
  //! Destructor
  virtual
  ~NoxObserver ();

  virtual
  void
  observeSolution(const Epetra_Vector& soln);
  
  void
  setStatisticsWriter( const Teuchos::RCP<Ginla::IO::StatsWriter>   & statsWriter,
                       const Teuchos::RCP<const Ginla::Operator::Virtual> & glOperator );
                       
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

    const ObserverType observerType_;
    const Teuchos::RCP<const Ginla::IO::StateWriter>          stateWriter_;
    const Teuchos::RCP<const Ginla::ModelEvaluator::Default>  modelEvaluator_;
    
    Teuchos::RCP<Ginla::IO::StatsWriter>   statsWriter_;
    Teuchos::RCP<const Ginla::Operator::Virtual> glOperator_;
};

} // namespace IO
} // namespace Ginla

#endif // GINLA_IO_NOXOBSERVER_H