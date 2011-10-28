/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

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

#include "VIO_Reader_Abstract.h"

#include <vtkCellArray.h>
#include <vtkPoints.h>

#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

#include <vtkFieldData.h>
#include <vtkPointData.h>

// =============================================================================
VIO::Reader::Abstract::
Abstract ( const std::string & filePath ) :
           filePath_ ( filePath )
{
}
// =============================================================================
VIO::Reader::Abstract::
~Abstract ()
{
}
// =============================================================================
double
VIO::Reader::Abstract::
extractDoubleValue_ ( const vtkSmartPointer<vtkDataArray> & array
                    ) const
{
    TEUCHOS_ASSERT_EQUALITY ( array->GetSize(), 1 );

    const vtkSmartPointer<vtkDoubleArray> arr =
        dynamic_cast<vtkDoubleArray*> ( array.GetPointer() );

    double val[1];
    vtkIdType i = 0;
    arr->GetTupleValue ( i, &*val );
    return val[0];
}
// =============================================================================
int
VIO::Reader::Abstract::
extractIntValue_ ( const vtkSmartPointer<vtkDataArray> & array
                 ) const
{
    TEUCHOS_ASSERT_EQUALITY ( array->GetSize(), 1 );

    const vtkSmartPointer<vtkIntArray> arr =
        dynamic_cast<vtkIntArray*> ( array.GetPointer() );

    int val[1];
    vtkIdType i = 0;
    arr->GetTupleValue ( i, &*val );
    return val[0];
}
// =============================================================================
Teuchos::Array<int>
VIO::Reader::Abstract::
extractIntArray_ ( const vtkSmartPointer<vtkDataArray> & array
                 ) const
{
    TEUCHOS_ASSERT_EQUALITY ( array->GetNumberOfComponents(), 1 );

    int n = array->GetNumberOfTuples();

    Teuchos::Array<int> intArray ( n );

    const vtkSmartPointer<vtkIntArray> arr =
        dynamic_cast<vtkIntArray*> ( array.GetPointer() );

    int val[1];
    for ( vtkIdType k=0; k<n; k++ )
    {
        arr->GetTupleValue ( k, &*val );
        intArray[k] = val[0];
    }

    return intArray;
}
// =============================================================================
Teuchos::Array<double>
VIO::Reader::Abstract::
extractDoubleArray_ ( const vtkSmartPointer<vtkDataArray> & array
                    ) const
{
    TEUCHOS_ASSERT_EQUALITY ( array->GetNumberOfComponents(), 1 );

    int n = array->GetNumberOfTuples();

    Teuchos::Array<double> doubleArray ( n );

    const vtkSmartPointer<vtkDoubleArray> arr =
        dynamic_cast<vtkDoubleArray*> ( array.GetPointer() );
    double val[1];
    for ( vtkIdType k=0; k<n; k++ )
    {
        arr->GetTupleValue ( k, &*val );
        doubleArray[k] = val[0];
    }

    return doubleArray;
}
// =============================================================================
Teuchos::ParameterList
VIO::Reader::Abstract::
readFieldData_ ( const vtkSmartPointer<vtkDataObject>  & dataObject
               ) const
{
    Teuchos::ParameterList fieldData;

    // extract the parameter values
    int numFieldData = dataObject->GetFieldData()->GetNumberOfArrays();
    for ( int k=0; k<numFieldData; k++ )
    {
        const vtkSmartPointer<vtkDataArray> & array =
            dataObject->GetFieldData()->GetArray ( k );

        std::string name = array->GetName();

        int numEntries = array->GetSize();

        int vtkDataTypeIndex = array->GetDataType();

        switch ( vtkDataTypeIndex )
        {
            // see http://www.vtk.org/doc/release/4.0/html/vtkSetGet_8h.html
        case 6:
            if ( numEntries==1 )
                fieldData.set ( name, this->extractIntValue_ ( array ) );
            else
                fieldData.set ( name, this->extractIntArray_ ( array ) );
            break;
        case 11:
            if ( numEntries==1 )
                fieldData.set ( name, this->extractDoubleValue_ ( array ) );
            else
                fieldData.set ( name, this->extractDoubleArray_ ( array ) );
            break;
        default:
            TEST_FOR_EXCEPTION ( true,
                                 std::runtime_error,
                                 "Illegal VTK data type index " <<
                                 vtkDataTypeIndex << "." );
        }
    }

    return fieldData;
}
// =============================================================================
