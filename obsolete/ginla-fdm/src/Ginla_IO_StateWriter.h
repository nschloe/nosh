/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl\"omer

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

#ifndef GINLA_IO_STATEWRITER_H
#define GINLA_IO_STATEWRITER_H

// Workaround for icpc's error "Include mpi.h before stdio.h"
#include <Teuchos_config.h>
#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include <LOCA_Parameter_Vector.H>

#include "Ginla_State_Virtual.h"

namespace Ginla {

namespace IO {

class StateWriter
{
public:
  //! Default constructor.
  StateWriter( const std::string & outputDir,
               const std::string & fileBaseName,
               const std::string & outputFormat,
               const unsigned int  maxIndex );

  void
  setOutputDir ( const string & directory );

  /** Writes a solution \c psi to a file with all parameters that
    * may be interesting.
    */
  void
  write ( const Teuchos::RCP<const Ginla::State::Virtual> & state,
          const unsigned int                                & index,
          const std::string                                 & filenameAppend,
          LOCA::ParameterVector                             & params
        ) const;

  void
  write ( const Teuchos::RCP<const Ginla::State::Virtual> & state,
          const unsigned int                                & index,
          LOCA::ParameterVector                             & params
        ) const;

  void
  write ( const Teuchos::RCP<const Ginla::State::Virtual> & state,
          const unsigned int                                & index,
          const std::string                                 & filenameAppend
        ) const;

  void
  write ( const Teuchos::RCP<const Ginla::State::Virtual> & state,
          const unsigned int                                & index
        ) const;

  protected:
  private:

  private:
    std::string outputDir_;
    std::string fileBaseName_;
    std::string filenameExtension_;
    unsigned int maxNumDigits_;
};

}

}

#endif // GINLA_IO_STATEWRITER_H
