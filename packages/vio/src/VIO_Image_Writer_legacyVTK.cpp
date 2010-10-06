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

#include "VIO_Image_Writer_legacyVTK.h"

#include <vtkStructuredPointsWriter.h>

// =============================================================================
// Constructor
VIO::Image::Writer::legacyVTK::
legacyVTK ( const std::string & filePath ) :
        VIO::Image::Writer::Abstract ( filePath ),
        isFormatBinary_( false ) // Set the default to ASCII.
                                 // If you want compressed files, use VTI.
{
}
// =============================================================================
// Destructor
VIO::Image::Writer::legacyVTK::
~legacyVTK()
{
}
// =============================================================================
void
VIO::Image::Writer::legacyVTK::
setFormatBinary()
{
  isFormatBinary_ = true;
}
// =============================================================================
void
VIO::Image::Writer::legacyVTK::
setFormatAscii()
{
  isFormatBinary_ = false;
}
// =============================================================================
void
VIO::Image::Writer::legacyVTK::
write () const
{
    // write the file
    vtkSmartPointer<vtkStructuredPointsWriter> writer =
        vtkSmartPointer<vtkStructuredPointsWriter>::New();
    writer->SetFileName ( filePath_.c_str() );
    writer->SetInput ( vtkDataSet_ );
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
