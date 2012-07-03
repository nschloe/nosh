// @HEADER
//
//    Helper class for writing out statistics.
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
#include "Cuantico_CsvWriter.hpp"

namespace Cuantico {
// ============================================================================
CsvWriter::
CsvWriter(const std::string &fileName,
          const std::string &delimeter
          ):
  fileStream_(),
  delimeter_( delimeter ),
  headerStart_("#"),
  doublePrec_(15),
  doubleColumnWidth_(doublePrec_ + 7),
  intColumnWidth_(5)
{
  // Set the output format
  // Think about replacing this with NOX::Utils::Sci.
  fileStream_.setf( std::ios::scientific );
  fileStream_.precision( 15 );

  fileStream_.open( fileName.c_str(), std::ios::trunc );
  return;
}
// ============================================================================
CsvWriter::
~CsvWriter()
{
  fileStream_.close();
}
// ============================================================================
void
CsvWriter::
writeHeader(const Teuchos::ParameterList & pList) const
{
  fileStream_ << headerStart_ << " ";

  bool isFirst = true;
  for (Teuchos::ParameterList::ConstIterator k=pList.begin(); k!=pList.end(); ++k )
  {
    if (!isFirst)
    {
      fileStream_ << delimeter_ << "  ";
      isFirst = false;
    }

    std::stringstream strstream;
    strstream.fill( ' ' );
    strstream << std::left;

    std::string label = pList.name( k );
    if ( pList.isType<double>( label ) )
      strstream.width( doubleColumnWidth_ );
    else if ( pList.isType<int>( label )
              || pList.isType<unsigned int>( label ) )
      strstream.width( (intColumnWidth_<label.length()) ?
                       label.length() : intColumnWidth_
                       );
    else if ( pList.isType<std::string>( label ) )
      strstream.width( pList.get<std::string>( label ).length() );
    else
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                           "Invalid data type for item \""
                           << label << "\"." );

    strstream << pList.name( k );

    // Write it out to the file.
    fileStream_ << strstream.str();
  }
  // flush:
  fileStream_ << std::endl;

  return;
}
// ============================================================================
void
CsvWriter::
writeRow(const Teuchos::ParameterList & pList) const
{
  // Pad from the left with the width of the header line starter
  // to get column alignment.
  fileStream_ << std::string(headerStart_.length() + 1, ' ');

  bool isFirst = true;
  for (Teuchos::ParameterList::ConstIterator k=pList.begin(); k!=pList.end(); ++k )
  {
    if (!isFirst)
    {
      fileStream_ << delimeter_ << "  ";
      isFirst = false;
    }

    std::stringstream strstream;
    strstream.fill( ' ' );
    strstream << std::left;    // have the number flush left

    std::string label = pList.name( k );
    if ( pList.isType<double>( label ) )
    {
      strstream.width( doubleColumnWidth_ );
      strstream.setf( std::ios::scientific );
      strstream.precision( doublePrec_ );
      strstream << pList.get<double>( label );
    }
    else if ( pList.isType<int>( label ) )
    {
      strstream.width( (intColumnWidth_<label.length()) ?
                       label.length() : intColumnWidth_
                       );
      strstream << pList.get<int>( label );
    }
    else if ( pList.isType<unsigned int>( label ) )
    {
      strstream.width( (intColumnWidth_<label.length()) ?
                       label.length() : intColumnWidth_
                       );
      strstream << pList.get<unsigned int>( label );
    }
    else if ( pList.isType<std::string>( label ) )
    {
      strstream.width( pList.get<std::string>( label ).length() );
      strstream << pList.get<std::string>( label );
    }
    else
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                           "Invalid data type for item \""
                           << label << "\"." );
  }

  // flush:
  fileStream_ << std::endl;

  return;
}
// ============================================================================
} // namespace Cuantico
