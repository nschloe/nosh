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
#include "Ginla_ModelEvaluator_Default.h"

// ============================================================================
Ginla::IO::NoxObserver::
NoxObserver ( const Teuchos::RCP<const Ginla::IO::StateWriter>         & stateWriter,
              const Teuchos::RCP<const Ginla::ModelEvaluator::Default> & modelEvaluator
            ) :
  stateWriter_ ( stateWriter ),
  modelEvaluator_ ( modelEvaluator )
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
observeSolution( const Epetra_Vector & soln )
{ 
    static int index = 0;
    index++;
    stateWriter_->write( modelEvaluator_->createState(soln),
                         index
                       );
    return;
}
// ============================================================================