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
extractDoubleValue ( const vtkSmartPointer<vtkDataArray> & array
                   ) const
{
    TEUCHOS_ASSERT_EQUALITY ( array->GetSize(), 1 );

    const vtkSmartPointer<vtkDoubleArray> arr = dynamic_cast<vtkDoubleArray*> ( array.GetPointer() );
    double val[1];
    vtkIdType i = 0;
    arr->GetTupleValue ( i, &*val );
    return val[0];
}
// =============================================================================
int
VIO::Reader::Abstract::
extractIntValue ( const vtkSmartPointer<vtkDataArray> & array
                ) const
{
    TEUCHOS_ASSERT_EQUALITY ( array->GetSize(), 1 );

    const vtkSmartPointer<vtkIntArray> arr = dynamic_cast<vtkIntArray*> ( array.GetPointer() );
    int val[1];
    vtkIdType i = 0;
    arr->GetTupleValue ( i, &*val );
    return val[0];
}
// =============================================================================
Teuchos::Array<int>
VIO::Reader::Abstract::
extractIntArray ( const vtkSmartPointer<vtkDataArray> & array
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
extractDoubleArray ( const vtkSmartPointer<vtkDataArray> & array
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
void
VIO::Reader::Abstract::
processImageData ( const vtkSmartPointer<vtkImageData>           & imageData,
                   Teuchos::RCP<ComplexMultiVector>              & z,
                   Teuchos::Array<int>                           & p,
                   UIntTuple                                     & dims,
                   DoubleTuple                                   & origin,
                   DoubleTuple                                   & spacing,
                   Teuchos::ParameterList                        & fieldData,
                   const Teuchos::RCP<const Teuchos::Comm<int> > & TComm
                 ) const
{
    double *tmp;
    tmp = imageData->GetOrigin();
    origin = Teuchos::tuple ( tmp[0], tmp[1] );
    tmp = imageData->GetSpacing();
    spacing = Teuchos::tuple ( tmp[0], tmp[1] );
    int *intTmp;
    intTmp = imageData->GetDimensions();
    dims = Teuchos::tuple ( ( unsigned int ) intTmp[0],
                            ( unsigned int ) intTmp[1] );

    fieldData = readFieldData ( imageData );

    // deep copy the permutation vector out, and delete the entry in the
    // parameter list
    p = fieldData.get<Teuchos::Array<int> > ( "p" );
    fieldData.remove ( "p" );

    vtkIdType numArrays = imageData->GetPointData()->GetNumberOfArrays();
    TEUCHOS_ASSERT_EQUALITY ( numArrays, 1 );

    const vtkSmartPointer<vtkDataArray> & array = imageData->GetPointData()->GetArray(0);
    
    vtkIdType numComponents = array->GetNumberOfComponents();
    TEUCHOS_ASSERT_EQUALITY ( numComponents, 2 );    // for *complex* values:

    // this is the total number of grid points, including the dummies
    // outside the domain
    vtkIdType numPoints = array->GetNumberOfTuples();

    // count the number of true grid points
    int numTruePoints = 0;
    for ( int k=0; k<numPoints; k++ )
        if ( p[k]>=0 )
            numTruePoints++;

    // create an appropriate map
    Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > ComplexMap =
        Teuchos::rcp ( new Tpetra::Map<Thyra::Ordinal> ( numTruePoints, 0, TComm ) );

    z = Teuchos::rcp ( new ComplexVector ( ComplexMap ) );

    // TODO only works on one core
    // fill z
    Teuchos::ArrayRCP<std::complex<double> > zView = z->get1dViewNonConst();
    
    double val[2];
    for ( int k = 0; k < numPoints; k++ )
    {
        if ( p[k]>=0 )
        {
            array->GetTuple(k,val);
            zView[p[k]] = std::complex<double> ( val[0], val[1] );
        }
    }

    return;
}
// =============================================================================
Teuchos::ParameterList
VIO::Reader::Abstract::
readFieldData ( const vtkSmartPointer<vtkImageData>  & imageData
              ) const
{
    Teuchos::ParameterList fieldData;

    // extract the parameter values
    int numFieldData = imageData->GetFieldData()->GetNumberOfArrays();
    for ( int k=0; k<numFieldData; k++ )
    {
        const vtkSmartPointer<vtkDataArray> & array =
            imageData->GetFieldData()->GetArray ( k );

        std::string name = array->GetName();

        int numEntries = array->GetSize();

        int vtkDataTypeIndex = array->GetDataType();

        switch ( vtkDataTypeIndex )
        {
            // see http://www.vtk.org/doc/release/4.0/html/vtkSetGet_8h.html
        case 6:
            if ( numEntries==1 )
                fieldData.set ( name, extractIntValue ( array ) );
            else
                fieldData.set ( name, extractIntArray ( array ) );
            break;
        case 11:
            if ( numEntries==1 )
                fieldData.set ( name, extractDoubleValue ( array ) );
            else
                fieldData.set ( name, extractDoubleArray ( array ) );
            break;
        default:
            TEST_FOR_EXCEPTION ( true,
                                 std::runtime_error,
                                 "Illegal VTK data type index " << vtkDataTypeIndex << "." );
        }
    }

    return fieldData;
}
// =============================================================================
