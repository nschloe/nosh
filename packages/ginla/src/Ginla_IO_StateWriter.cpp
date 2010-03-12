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

#include "Ginla_IO_StateWriter.h"

#include "Ginla_Helpers.h"

// ============================================================================
Ginla::IO::StateWriter::
StateWriter( const std::string & outputDir,
             const std::string & fileBaseName,
             const std::string & outputFormat,
             const unsigned int maxIndex ):
    outputDir_( outputDir ),
    fileBaseName_( fileBaseName )
{
  if ( outputFormat.compare("VTI")==0 )
    filenameExtension_ = "vti";
  else if ( outputFormat.compare("VTK")==0 )
    filenameExtension_ = "vtk";
  else
    TEST_FOR_EXCEPTION( true,
                        std::runtime_error,
                        "outputFormat_ (\"" << outputFormat
                        << "\") must be either one of \"VTI\", \"VTK\"" );

  maxNumDigits_ = numDigits( maxIndex );
  return;
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
write ( const Teuchos::RCP<const ComplexVector>        & psi,
        const Teuchos::RCP<const Recti::Grid::General> & grid,
        const unsigned int                             & index,
        const std::string                              & filenameAppend,
        LOCA::ParameterVector                          & params
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
    grid->writeWithGrid ( *psi, *p, fileName.str() );

    return;
}
// ============================================================================
void
Ginla::IO::StateWriter::
write ( const Teuchos::RCP<const ComplexVector>        & psi,
        const Teuchos::RCP<const Recti::Grid::General> & grid,
        const unsigned int                             & index,
        const std::string                              & filenameAppend
      ) const
{
  LOCA::ParameterVector empty;
  write ( psi, grid, index, filenameAppend, empty );
  return;
}
// ============================================================================
void
Ginla::IO::StateWriter::
write ( const Teuchos::RCP<const ComplexVector>        & psi,
        const Teuchos::RCP<const Recti::Grid::General> & grid,
        const unsigned int                             & index,
        LOCA::ParameterVector                          & params
      ) const
{
  std::string filenameAppend = "";
  write ( psi, grid, index, filenameAppend, params );
  return;
}
// ============================================================================
void
Ginla::IO::StateWriter::
write ( const Teuchos::RCP<const ComplexVector>        & psi,
        const Teuchos::RCP<const Recti::Grid::General> & grid,
        const unsigned int                             & index
      ) const
{
  std::string filenameAppend = "";
  LOCA::ParameterVector empty;
  write ( psi, grid, index, filenameAppend, empty );
  return;
}
// ============================================================================
unsigned int
Ginla::IO::StateWriter::
numDigits ( const int i )
{
    int numDigits = 0;
    int ii = i;
    if ( ii < 0 )
        ii = -ii;

    while ( ii > 0 )
    {
        numDigits++;
        ii/=10;
    }
    return numDigits;
}
// ============================================================================