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

#include "VtkWriter.h"

#include <vtkStructuredPointsWriter.h>

// =============================================================================
// Constructor
VtkWriter::VtkWriter ( const std::string & filePath ) :
        AbstractImageWriter ( filePath ),
        isFormatBinary_( false ) // Set the default to ASCII.
                                 // If you want compressed files, use VTI.
{
}
// =============================================================================
// Destructor
VtkWriter::~VtkWriter()
{
}
// =============================================================================
void
VtkWriter::setFormatBinary()
{
  isFormatBinary_ = true;
}
// =============================================================================
void
VtkWriter::setFormatAscii()
{
  isFormatBinary_ = false;
}
// =============================================================================
void
VtkWriter::write () const
{
    // write the file
    vtkSmartPointer<vtkStructuredPointsWriter> writer =
        vtkSmartPointer<vtkStructuredPointsWriter>::New();
    writer->SetFileName ( filePath_.c_str() );
    writer->SetInput ( imageData_ );
    if ( isFormatBinary_ )
        // write binary data to avoid losing precision
        writer->SetFileTypeToBinary();
    else
        // ASCII format stores only eight significant digits
        writer->SetFileTypeToASCII();

    writer->Write();

    return;
}
// =============================================================================
