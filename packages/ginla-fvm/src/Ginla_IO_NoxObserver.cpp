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

#include "Ginla_IO_NoxObserver.h"

#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_StatsWriter.h"
#include "Ginla_Helpers.h"
#include "Ginla_EpetraFVM_ModelEvaluator.h"

// ============================================================================
Ginla::IO::NoxObserver::
NoxObserver ( const Teuchos::RCP<const Ginla::EpetraFVM::ModelEvaluator> & modelEval,
              const NoxObserver::ObserverType                            & observerType
            ) :
  modelEval_ ( modelEval ),
  observerType_( observerType ),
  statsWriter_ ( Teuchos::null )
{
}
// ============================================================================
Ginla::IO::NoxObserver::
~NoxObserver ()
{
}
// ============================================================================
void
Ginla::IO::NoxObserver::
setStatisticsWriter( const Teuchos::RCP<Ginla::IO::StatsWriter> & statsWriter
                   )
{
  statsWriter_   = statsWriter;
  return;
}
// ============================================================================
void
Ginla::IO::NoxObserver::
observeSolution( const Epetra_Vector & soln )
{
    // define state
    const Teuchos::RCP<const Ginla::EpetraFVM::State> savable =
        modelEval_->createSavable( soln );

    // The switch hack is necessary as different continuation algorithms
    // call printSolution() a different number of times per step, e.g.,
    // to store solutions, null vectors, and so forth.
    switch ( observerType_ )
    {
      case NONLINEAR:
          savable->save( 0 );
          break;
      case CONTINUATION:
          this->observeContinuation_( savable );
          break;
      case TURNING_POINT:
          this->observeTurningPointContinuation_( savable );
          break;
      default:
          TEST_FOR_EXCEPTION ( true,
                               std::logic_error,
                               "Illegal observer type " << observerType_ );
    }

    return;
}
// ============================================================================
void
Ginla::IO::NoxObserver::
observeContinuation_( const Teuchos::RCP<const Ginla::EpetraFVM::State> & state )
{
  static int index = -1;
  index++;

  this->saveContinuationStatistics_( index, state );

  TEUCHOS_ASSERT( !modelEval_.is_null() );
  Teuchos::RCP<LOCA::ParameterVector> p = modelEval_->getParameters();

  Teuchos::RCP<Teuchos::ParameterList> tp = Ginla::Helpers::locaParameterVector2teuchosParameterList( *p );
  state->save( index, *tp );

  return;
}
// ============================================================================
void
Ginla::IO::NoxObserver::
observeTurningPointContinuation_( const Teuchos::RCP<const Ginla::EpetraFVM::State> & state )
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
        TEST_FOR_EXCEPTION( true,
                            std::logic_error,
                            "Not yet implemented." <<
                            "This part of the code used to write state and null vector alternately for turning point continuation" <<
                            "but because of how StkMesh is organized, it seems impossible to first write to one file, then to" <<
                            "another with with the same mesh. Need to investigate." );
    return;
}
// ============================================================================
void
Ginla::IO::NoxObserver::
saveContinuationStatistics_( const int stepIndex,
                             const Teuchos::RCP<const Ginla::EpetraFVM::State> & state
                           )
{
    if ( !statsWriter_.is_null() )
    {
        TEUCHOS_ASSERT( !state.is_null() );
        Teuchos::RCP<Teuchos::ParameterList> paramList = statsWriter_->getListNonConst();

        paramList->set( "0step", stepIndex );

        // put the parameter list into statsWriter_
        std::string labelPrepend = "1";
        TEUCHOS_ASSERT( !modelEval_.is_null() );
        Ginla::Helpers::appendToTeuchosParameterList( *paramList,
                                                      *(modelEval_->getParameters()),
                                                      labelPrepend
                                                    );

        paramList->set( "2free energy", state->freeEnergy() );
        paramList->set( "2||x||_2 scaled", state->normalizedScaledL2Norm() );
        paramList->set( "2vorticity", state->getVorticity() );

        // actually print the data
        statsWriter_->print();
    }
}
// ============================================================================
