/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"oemer

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

#include "VIO_Writer_Abstract.h"

#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkFieldData.h>
#include <vtkDoubleArray.h>

// =============================================================================
VIO::Writer::Abstract::
Abstract ( const std::string & filePath ) :
        filePath_ ( filePath ),
        vtkDataSet_( 0 )
{
}
// =============================================================================
VIO::Writer::Abstract::
~Abstract()
{
}
// =============================================================================
void
VIO::Writer::Abstract::
addFieldData ( const Teuchos::Array<int> & array,
               const std::string         & name
             )
{
    // create field data
    vtkSmartPointer<vtkIntArray> fieldData = vtkSmartPointer<vtkIntArray>::New();

    fieldData->SetName ( name.c_str() );

    // fill the field
    for ( int k=0; k<array.length(); k++ )
        fieldData->InsertNextValue ( array[k] );

    vtkDataSet_->GetFieldData()->AddArray ( fieldData );

    return;
}
// =============================================================================
void
VIO::Writer::Abstract::
addParameterList ( const Teuchos::ParameterList & problemParams )
{
    // add to imageData_
    Teuchos::map<std::string, Teuchos::ParameterEntry>::const_iterator k;
    for ( k = problemParams.begin(); k != problemParams.end(); ++k )
    {
        std::string paramName = problemParams.name ( k );
        if ( problemParams.isType<int> ( paramName ) )
        {
            vtkSmartPointer<vtkIntArray> fieldData = vtkSmartPointer<vtkIntArray>::New();
            fieldData->InsertNextValue ( problemParams.get<int> ( paramName ) );
            fieldData->SetName ( paramName.c_str() );
            vtkDataSet_->GetFieldData()->AddArray ( fieldData );
        }
        else if ( problemParams.isType<double> ( paramName ) )
        {
            vtkSmartPointer<vtkDoubleArray> fieldData = vtkSmartPointer<vtkDoubleArray>::New();
            fieldData->InsertNextValue ( problemParams.get<double> ( paramName ) );
            fieldData->SetName ( paramName.c_str() );
            vtkDataSet_->GetFieldData()->AddArray ( fieldData );
        }
        else if ( problemParams.isType<Teuchos::Array<double> > ( paramName ) )
        {
            vtkSmartPointer<vtkDoubleArray> fieldData = vtkSmartPointer<vtkDoubleArray>::New();
            const Teuchos::Array<double> & arr = problemParams.get<Teuchos::Array<double> > ( paramName );
            for ( int k=0; k<arr.length(); k++ )
                fieldData->InsertNextValue ( arr[k] );
            fieldData->SetName ( paramName.c_str() );
            vtkDataSet_->GetFieldData()->AddArray ( fieldData );
        }
        else
        {
            TEST_FOR_EXCEPTION ( true,
                                 std::runtime_error,
                                 "Illegal type of parameter \"" << paramName << "\"." );
        }
    }

    return;
}
// =============================================================================