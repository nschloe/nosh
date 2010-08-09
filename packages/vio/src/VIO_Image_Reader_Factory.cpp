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

#include "VIO_Image_Reader_Factory.h"

#include "VIO_Image_Reader_VTI.h"
#include "VIO_Image_Reader_legacyVTK.h"

// =============================================================================
Teuchos::RCP<VIO::Image::Reader::Abstract>
VIO::Image::Reader::Factory::
create ( const std::string & fileName )
{
    // analyze the file name for extension
    int         dotPos    = fileName.rfind ( "." );
    std::string extension = fileName.substr ( dotPos+1, fileName.size()-dotPos-1 );

    if ( extension.compare ( "vtk" ) == 0 )
    {
        return Teuchos::rcp( new VIO::Image::Reader::legacyVTK ( fileName ) );
    }
    else if ( extension.compare ( "vti" ) == 0 )
    {
        return Teuchos::rcp( new VIO::Image::Reader::VTI ( fileName ) );
    }
    else
    {
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Error when reading file \"" << fileName
                             << "\". File name extension \"" << extension << "\" "
                             << "not recognized. Must be one of \"vtk\", "
                             << "\"vti\"." );
    }
}
// =============================================================================