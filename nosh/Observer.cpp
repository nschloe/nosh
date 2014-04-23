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

#include "nosh/Observer.hpp"

#include <string>

#include "nosh/ModelEvaluator_Virtual.hpp"
#include "nosh/StkMesh.hpp"

namespace Nosh
{
// ============================================================================
Observer::
Observer(const Teuchos::RCP<const Nosh::ModelEvaluator::Virtual> &modelEval,
         const std::string & csvFilename,
         const std::string & contParamName,
         const bool isTurningPointContinuation
        ) :
  modelEval_(modelEval),
  csvWriter_(csvFilename, " "),
  contParamName_(contParamName),
  isTurningPointContinuation_(isTurningPointContinuation)
{
}
// ============================================================================
Observer::
~Observer ()
{
}
// ============================================================================
void
Observer::
observeSolution(const Epetra_Vector &soln)
{
  modelEval_->getMesh()->write(soln, 0.0);

  return;
}
// ============================================================================
void
Observer::
observeSolution(const Epetra_Vector& soln,
                double paramVal)
{
  // This if-else hack is necessary as different continuation algorithms
  // call printSolution() a different number of times per step, e.g.,
  // to store solutions, null vectors, and so forth.
  if (isTurningPointContinuation_)
    this->observeTurningPointContinuation_(soln, paramVal);
  else
    this->observeContinuation_(soln, paramVal);

  return;
}
// ============================================================================
void
Observer::
observeContinuation_(const Epetra_Vector &soln,
                     const double paramVal
                    )
{
  static int index = -1;
  index++;

  this->saveContinuationStatistics_(soln, paramVal, index);

  // Storing the parameter value as "time" variable here is convenient, but
  // has a downside: The default output format ExodusII insists that the
  // values for time are monotonically increasing. The parameter, however,
  // can decrease. Since there's no hard reason for the monotonicity condition,
  // many things will continue to work fine if the time data isn't monotonous.
  // The display in ParaView is one example where it doesn't work so well.
  // As a work-around for that, paramVal could be replaced by index.
  modelEval_->getMesh()->write(soln, paramVal);

  return;
}
// ============================================================================
void
Observer::
observeTurningPointContinuation_(const Epetra_Vector &soln,
                                 const double paramVal
                                )
{
  static int index = -1;
  static bool isSolution = false;

  // alternate between solution and nullvector
  isSolution = !isSolution;
  if ( isSolution ) {
    index++;
    this->saveContinuationStatistics_(soln, paramVal, index);
    modelEval_->getMesh()->write(soln, index);
  } else
    TEUCHOS_TEST_FOR_EXCEPT_MSG( true, "Not yet implemented." );
  // This part of the code used to write state and null vector alternately
  // for turning point continuation, but because of how StkMesh is
  // organized, it seems impossible to first write to one file, then to
  // another with with the same mesh. Need to investigate.
  return;
}
// ============================================================================
void
Observer::
saveContinuationStatistics_(const Epetra_Vector &soln,
                            const double paramVal,
                            const int stepIndex
                           )
{
  // Construct parameter list to stuff into the csvWriter_.
  Teuchos::ParameterList paramList;
  paramList.set( "(0) step", stepIndex );

  // Continuation parameter.
  paramList.set("(1) "+contParamName_, paramVal);

  // Some extra stats.
  paramList.set( "(2) Gibbs energy", modelEval_->gibbsEnergy(soln) );
  paramList.set( "(2) ||x||_2 scaled", modelEval_->norm(soln) );

  // Write out header.
  if (stepIndex == 0)
    csvWriter_.writeHeader(paramList);
  // Write out the data.
  csvWriter_.writeRow(paramList);
  return;
}
// ============================================================================
} // namespace Nosh
