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

#include "VIO_Image_Writer_VTI.h"

#include <vtkXMLImageDataWriter.h>

// =============================================================================
// Constructor
VIO::Image::Writer::VTI::
VTI ( const std::string & filePath ) :
      VIO::Image::Writer::Abstract ( filePath )
{
}
// =============================================================================
// Destructor
VIO::Image::Writer::VTI::
~VTI()
{
}
// =============================================================================
void
VIO::Image::Writer::VTI::
write () const
{
    // write the file
    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName ( filePath_.c_str() );
    writer->SetInput ( imageData_ );
    writer->Write();

    return;
}
// =============================================================================
