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

#include "Ginla_NoxObserver.hpp"

#include "Ginla_State.hpp"
#include "Ginla_StatsWriter.hpp"
#include "Ginla_Helpers.hpp"
#include "Ginla_ModelEvaluator.hpp"

namespace Ginla {
// ============================================================================
NoxObserver::
NoxObserver ( const Teuchos::RCP<const Ginla::ModelEvaluator> & modelEval,
              const NoxObserver::EObserverType                & observerType
            ) :
  modelEval_ ( modelEval ),
  observerType_( observerType ),
  statsWriter_ ( Teuchos::null )
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
setStatisticsWriter( const Teuchos::RCP<Ginla::StatsWriter> & statsWriter
                   )
{
  statsWriter_ = statsWriter;
  return;
}
// ============================================================================
void
NoxObserver::
observeSolution( const Epetra_Vector & soln )
{
    // define state
    const Teuchos::RCP<const Ginla::State> savable =
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
          TEST_FOR_EXCEPT_MSG( true,
                                       "Illegal observer type " << observerType_ );
    }

    return;
}
// ============================================================================
void
NoxObserver::
observeContinuation_( const Teuchos::RCP<const Ginla::State> & state )
{
  static int index = -1;
  index++;

  this->saveContinuationStatistics_( index, state );

#ifdef _DEBUG_
  TEUCHOS_ASSERT( !state.is_null() );
#endif
  state->save( index );

  return;
}
// ============================================================================
void
NoxObserver::
observeTurningPointContinuation_( const Teuchos::RCP<const Ginla::State> & state )
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
        TEST_FOR_EXCEPT_MSG( true, "Not yet implemented." );
        // This part of the code used to write state and null vector alternately
        // for turning point continuation, but because of how StkMesh is
        // organized, it seems impossible to first write to one file, then to
        // another with with the same mesh. Need to investigate.
    return;
}
// ============================================================================
void
NoxObserver::
saveContinuationStatistics_( const int stepIndex,
                             const Teuchos::RCP<const Ginla::State> & state
                           )
{
    if ( !statsWriter_.is_null() )
    {
#ifdef _DEBUG_
        TEUCHOS_ASSERT( !state.is_null() );
#endif
        Teuchos::RCP<Teuchos::ParameterList> paramList =
            statsWriter_->getListNonConst();

        paramList->set( "0step", stepIndex );

        // put the parameter list into statsWriter_
        std::string labelPrepend = "1";
#ifdef _DEBUG_
        TEUCHOS_ASSERT( !modelEval_.is_null() );
#endif
        Ginla::Helpers::appendToTeuchosParameterList( *paramList,
                                                      *(modelEval_->getParameters()),
                                                      labelPrepend
                                                    );

        paramList->set( "2free energy", state->freeEnergy() );
        paramList->set( "2||x||_2 scaled", state->normalizedScaledL2Norm() );

        // actually print the data
        statsWriter_->print();
    }
}
// ============================================================================
} // namespace Ginla
