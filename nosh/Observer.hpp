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
#include <string>

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <NOX_Epetra_Observer.H>

#include "nosh/CsvWriter.hpp"
// =============================================================================
// forward declarations
namespace Komplex2
{
class LinearProblem;
}
namespace Nosh
{
namespace ModelEvaluator
{
class Virtual;
}
}
// =============================================================================
namespace Nosh
{

class Observer:  public NOX::Epetra::Observer
{
public:
//! Constructor
  Observer(const std::shared_ptr<const Nosh::ModelEvaluator::Virtual> &modelEval,
           const std::string & csvFilename = "",
           const std::string & contParamName = "",
           const bool isTurningPointContinuation = false
          );

//! Destructor
  virtual
  ~Observer ();

  virtual
  void
  observeSolution(const Epetra_Vector &soln);

  virtual
  void
  observeSolution(const Epetra_Vector& soln,
                  double paramVal);

protected:
private:
  void
  observeContinuation_(const Epetra_Vector &soln,
                       const double paramVal
                     );

  void
  observeTurningPointContinuation_(const Epetra_Vector &soln,
                                   const double paramVal
                                 );

  void
  saveContinuationStatistics_(const Epetra_Vector &soln,
                              const double paramVal,
                              const int stepIndex
                             );

private:
  const std::shared_ptr<const Nosh::ModelEvaluator::Virtual> modelEval_;
  Nosh::CsvWriter csvWriter_;
  const std::string contParamName_;
  const bool isTurningPointContinuation_;
};
} // namespace Nosh
#endif // NOSH_NOXOBSERVER_H
