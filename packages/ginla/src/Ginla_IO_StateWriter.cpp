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

#include "Ginla_IO_StateWriter.h"

#include "Ginla_Helpers.h"

// ============================================================================
Ginla::IO::StateWriter::
StateWriter( const std::string & outputDir,
             const std::string & fileBaseName,
             const std::string & filenameExtension,
             const unsigned int maxIndex ):
    outputDir_( outputDir ),
    fileBaseName_( fileBaseName ),
    filenameExtension_( filenameExtension ),
    maxNumDigits_( Ginla::Helpers::numDigits( maxIndex ) )
{
}
// ============================================================================
void
Ginla::IO::StateWriter::
setOutputDir ( const string & directory )
{
  outputDir_ = directory;
  return;
}
// ============================================================================
void
Ginla::IO::StateWriter::
write ( const Teuchos::RCP<const Ginla::State::Virtual> & state,
        const unsigned int                              & index,
        const std::string                               & filenameAppend,
        LOCA::ParameterVector                           & params
      ) const
{
    // get the parameter list
    Teuchos::RCP<Teuchos::ParameterList> p =
        Ginla::Helpers::locaParameterVector2teuchosParameterList( params );        

    // create the file name
    stringstream fileName;
    fileName
    << outputDir_ << "/" << fileBaseName_
    << setw ( maxNumDigits_ ) << setfill ( '0' ) << index
    << filenameAppend << "." << filenameExtension_;

    // write the file
    state->save( fileName.str(), *p );

    return;
}
// ============================================================================
void
Ginla::IO::StateWriter::
write ( const Teuchos::RCP<const Ginla::State::Virtual> & state,
        const unsigned int                              & index,
        const std::string                               & filenameAppend
      ) const
{
  LOCA::ParameterVector empty;
  this->write ( state, index, filenameAppend, empty );
  return;
}
// ============================================================================
void
Ginla::IO::StateWriter::
write ( const Teuchos::RCP<const Ginla::State::Virtual> & state,
        const unsigned int                              & index,
        LOCA::ParameterVector                           & params
      ) const
{
  std::string filenameAppend = "";
  this->write ( state, index, filenameAppend, params );
  return;
}
// ============================================================================
void
Ginla::IO::StateWriter::
write ( const Teuchos::RCP<const Ginla::State::Virtual> & state,
        const unsigned int                              & index
      ) const
{
  std::string filenameAppend = "";
  LOCA::ParameterVector empty;
  this->write ( state, index, filenameAppend, empty );
  return;
}
// ============================================================================