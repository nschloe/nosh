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

#include "Cuantico_NoxObserver.hpp"

#include "Cuantico_State.hpp"
#include "Cuantico_CsvWriter.hpp"
#include "Cuantico_Helpers.hpp"
#include "Cuantico_ModelEvaluator.hpp"

namespace Cuantico {
// ============================================================================
NoxObserver::
NoxObserver(const Teuchos::RCP<const Cuantico::ModelEvaluator> &modelEval,
            const std::string & filename,
            const NoxObserver::EObserverType &observerType
             ) :
  modelEval_(modelEval),
  csvWriter_(new Cuantico::CsvWriter(filename, " ")),
  observerType_(observerType)
{
}
// ============================================================================
NoxObserver::
~NoxObserver ()
{
}
// ============================================================================
void
NoxObserver::
observeSolution(const Epetra_Vector &soln)
{
  // define state
  const Teuchos::RCP<const Cuantico::State> savable =
    modelEval_->createSavable( soln );

  // The switch hack is necessary as different continuation algorithms
  // call printSolution() a different number of times per step, e.g.,
  // to store solutions, null vectors, and so forth.
  switch ( observerType_ )
  {
    case OBSERVER_TYPE_NEWTON:
      savable->save( 0 );
      break;
    case OBSERVER_TYPE_CONTINUATION:
      this->observeContinuation_( savable );
      break;
    case OBSERVER_TYPE_TURNING_POINT:
      this->observeTurningPointContinuation_( savable );
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                           "Illegal observer type " << observerType_ );
  }

  return;
}
// ============================================================================
void
NoxObserver::
observeContinuation_( const Teuchos::RCP<const Cuantico::State> &state )
{
  static int index = -1;
  index++;

  this->saveContinuationStatistics_(index, state);

#ifdef _DEBUG_
  TEUCHOS_ASSERT( !state.is_null() );
#endif
  state->save( index );

  return;
}
// ============================================================================
void
NoxObserver::
observeTurningPointContinuation_( const Teuchos::RCP<const Cuantico::State> &state )
{
  static int index = -1;
  static bool isSolution = false;

  // alternate between solution and nullvector
  isSolution = !isSolution;
  if ( isSolution )
  {
    index++;
    this->saveContinuationStatistics_( index, state );
    state->save( index );
  }
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG( true, "Not yet implemented." );
  // This part of the code used to write state and null vector alternately
  // for turning point continuation, but because of how StkMesh is
  // organized, it seems impossible to first write to one file, then to
  // another with with the same mesh. Need to investigate.
  return;
}
// ============================================================================
void
NoxObserver::
saveContinuationStatistics_(const int stepIndex,
                            const Teuchos::RCP<const Cuantico::State> &state
                            )
{
  if ( !csvWriter_.is_null() )
  {
#ifdef _DEBUG_
    TEUCHOS_ASSERT( !state.is_null() );
#endif
    // Construct parameter list to stuff into the csvWriter_.
    Teuchos::ParameterList paramList;
    paramList.set( "0step", stepIndex );
    // Model evaluator parameters.
    std::string labelPrepend = "1";
#ifdef _DEBUG_
    TEUCHOS_ASSERT( !modelEval_.is_null() );
    TEUCHOS_ASSERT( !modelEval_->get_p_latest().is_null() );
#endif
    Teuchos::RCP<const Epetra_Vector> meParams =
      modelEval_->get_p_latest();
    // TODO cache this
    Teuchos::RCP<const Teuchos::Array<std::string> > names =
      modelEval_->get_p_names(0);
    for (int k=0; k<names->length(); k++)
      paramList.set("1"+(*names)[k], (*meParams)[k]);
    // Some extra stats.
    paramList.set( "2free energy", state->freeEnergy() );
    paramList.set( "2||x||_2 scaled", state->normalizedScaledL2Norm() );

    // Write out header.
    if (stepIndex == 0)
      csvWriter_->writeHeader(paramList);
    // Write out the data.
    csvWriter_->writeRow(paramList);
  }
  return;
}
// ============================================================================
} // namespace Cuantico
