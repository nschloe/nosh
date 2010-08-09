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

#include "VIO_Image_Reader_legacyVTK.h"

#include <vtkStructuredPointsReader.h>
#include <vtkPointData.h>

// =============================================================================
VIO::Image::Reader::legacyVTK::
legacyVTK ( const std::string & filePath ) :
        VIO::Image::Reader::Abstract ( filePath )
{
}
// =============================================================================
VIO::Image::Reader::legacyVTK::
~legacyVTK ()
{
}
// =============================================================================
void
VIO::Image::Reader::legacyVTK::
read ( Teuchos::RCP<ComplexMultiVector>              & z,
       Teuchos::Array<int>                           & p,
       UIntTuple                                     & dims,
       DoubleTuple                                   & origin,
       DoubleTuple                                   & spacing,
       Teuchos::ParameterList                        & fieldData,
       const Teuchos::RCP<const Teuchos::Comm<int> > & TComm
     ) const
{
    vtkSmartPointer<vtkStructuredPointsReader> reader =
        vtkSmartPointer<vtkStructuredPointsReader>::New();
    reader->SetFileName ( filePath_.c_str() );
    reader->Update();
    
    processImageData ( reader->GetOutput(),
                       z, p, dims, origin, spacing, fieldData, TComm );

    return;
}
// =============================================================================