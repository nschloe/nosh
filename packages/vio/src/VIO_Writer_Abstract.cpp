/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distri7buted in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "VIO_Writer_Abstract.h"

#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <EpetraExt_Utils.h> // for toString

// =============================================================================
VIO::Writer::Abstract::
Abstract ( const std::string & filePath ) :
        filePath_ ( filePath ),
        imageData_ ( vtkSmartPointer<vtkImageData>::New() )
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
               const std::string         & name )
{
    // create field data
    vtkSmartPointer<vtkIntArray> fieldData = vtkSmartPointer<vtkIntArray>::New();

    fieldData->SetName ( name.c_str() );

    // fill the field
    for ( int k=0; k<array.length(); k++ )
        fieldData->InsertNextValue ( array[k] );

    imageData_->GetFieldData()->AddArray ( fieldData );

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
            imageData_->GetFieldData()->AddArray ( fieldData );
        }
        else if ( problemParams.isType<double> ( paramName ) )
        {
            vtkSmartPointer<vtkDoubleArray> fieldData = vtkSmartPointer<vtkDoubleArray>::New();
            fieldData->InsertNextValue ( problemParams.get<double> ( paramName ) );
            fieldData->SetName ( paramName.c_str() );
            imageData_->GetFieldData()->AddArray ( fieldData );
        }
        else if ( problemParams.isType<Teuchos::Array<double> > ( paramName ) )
        {
            vtkSmartPointer<vtkDoubleArray> fieldData = vtkSmartPointer<vtkDoubleArray>::New();
            const Teuchos::Array<double> & arr = problemParams.get<Teuchos::Array<double> > ( paramName );
            for ( int k=0; k<arr.length(); k++ )
                fieldData->InsertNextValue ( arr[k] );
            fieldData->SetName ( paramName.c_str() );
            imageData_->GetFieldData()->AddArray ( fieldData );
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
void
VIO::Writer::Abstract::
setImageData ( const Epetra_MultiVector              & x,
               const Teuchos::Tuple<unsigned int,2>  & Nx,
               const Teuchos::Tuple<double,2>        & h,
               const Teuchos::Array<int>             & p,
               const Teuchos::Array<std::string>     & scalarsNames
             )
{
    int numVecs   = x.NumVectors();
    int numPoints = ( Nx[0]+1 ) * ( Nx[1]+1 );

    // get scalarsNames, and insert default names if empty
    Teuchos::Array<std::string> scNames ( scalarsNames );
    if ( scNames.empty() )
    {
        scNames.resize ( numVecs );
        for ( int vec=0; vec<numVecs; vec++ )
            scNames[vec] = "x" + EpetraExt::toString ( vec );
    }

    // set other image data
    imageData_->SetDimensions ( Nx[0]+1, Nx[1]+1, 1 );
    imageData_->SetOrigin ( 0.0, 0.0, 0.0 );
    imageData_->SetSpacing ( h[0], h[1], 0.0 );

    // fill the scalar field
    vtkSmartPointer<vtkDoubleArray> scalars =
        vtkSmartPointer<vtkDoubleArray>::New();

    bool isScrambled = !p.empty();

    if ( isScrambled )
    {
        TEUCHOS_ASSERT_EQUALITY ( numPoints, p.length() );
        addFieldData ( p, "p" );
    }

    // fill the scalars vector and add it to imageData_
    if ( isScrambled )
    {
        double dummy = 0.0;
        for ( int vec=0; vec<numVecs; vec++ )
        {
            scalars->SetName ( scNames[vec].c_str() );
            for ( int k=0; k<numPoints; k++ )
                scalars->InsertNextValue ( p[k]>=0 ? x[vec][p[k]] : dummy );
            imageData_->GetPointData()->AddArray ( scalars );
        }
    }
    else
        for ( int vec=0; vec<numVecs; vec++ )
        {
            scalars->SetName ( scNames[vec].c_str() );
            for ( int k=0; k<numPoints; k++ )
                scalars->InsertNextValue ( x[vec][k] );
            imageData_->GetPointData()->AddArray ( scalars );
        }

    return;
}
// =============================================================================
void
VIO::Writer::Abstract::
setImageData ( const DoubleMultiVector               & x,
               const Teuchos::Tuple<unsigned int,2>  & Nx,
               const Teuchos::Tuple<double,2>        & h,
               const Teuchos::Array<int>             & p,
               const Teuchos::Array<std::string>     & scalarsNames
             )
{
    int numVecs   = x.getNumVectors();
    int numPoints = ( Nx[0]+1 ) * ( Nx[1]+1 );

    // get scalarsNames, and insert default names if empty
    Teuchos::Array<std::string> scNames ( scalarsNames );
    if ( scNames.empty() )
    {
        scNames.resize ( numVecs );
        for ( int vec=0; vec<numVecs; vec++ )
            scNames[vec] = "x" + EpetraExt::toString ( vec );
    }

    // set other image data
    imageData_->SetDimensions ( Nx[0]+1, Nx[1]+1, 1 );
    imageData_->SetOrigin ( 0.0, 0.0, 0.0 );
    imageData_->SetSpacing ( h[0], h[1], 0.0 );

    // fill the scalar field
    vtkSmartPointer<vtkDoubleArray> scalars =
        vtkSmartPointer<vtkDoubleArray>::New();

    double dummy = 0.0;
    bool isScrambled = !p.empty();
    if ( isScrambled )
    {
        TEUCHOS_ASSERT_EQUALITY ( numPoints, p.length() );
        addFieldData ( p, "p" );
    }

    // fill the scalars vector and add it to imageData_
    Teuchos::ArrayRCP<const double> xView;
    for ( int vec=0; vec<numVecs; vec++ )
    {
        xView = x.getVector ( vec )->get1dView();
        scalars->SetName ( scNames[vec].c_str() );
        for ( int k=0; k<numPoints; k++ )
        {
            if ( isScrambled )
                scalars->InsertNextValue ( p[k]>=0 ? xView[p[k]] : dummy );
            else
                scalars->InsertNextValue ( xView[k] );
        }
        imageData_->GetPointData()->AddArray ( scalars );
    }

    return;
}
// =============================================================================
void
VIO::Writer::Abstract::
setImageData ( const ComplexMultiVector              & x,
               const Teuchos::Tuple<unsigned int,2>  & Nx,
               const Teuchos::Tuple<double,2>        & h,
               const Teuchos::Array<int>             & p,
               const Teuchos::Array<std::string>     & scalarsNames
             )
{
    int numVecs   = x.getNumVectors();
    int numPoints = ( Nx[0]+1 ) * ( Nx[1]+1 );

    // get scalarsNames, and insert default names if empty
    Teuchos::Array<std::string> scNames ( scalarsNames );
    if ( scNames.empty() )
    {
        scNames.resize ( numVecs );
        for ( int vec=0; vec<numVecs; vec++ )
            scNames[vec] = "z" + EpetraExt::toString ( vec );
    }

    // set other image data
    imageData_->SetDimensions ( Nx[0]+1, Nx[1]+1, 1 );
    imageData_->SetOrigin ( 0.0, 0.0, 0.0 );
    imageData_->SetSpacing ( h[0], h[1], 0.0 );

    // fill the scalar field
    vtkSmartPointer<vtkDoubleArray> scalars =
        vtkSmartPointer<vtkDoubleArray>::New();

    double dummy = 0.0;
    bool isScrambled = !p.empty();
    if ( isScrambled )
    {
        TEUCHOS_ASSERT_EQUALITY ( numPoints, p.length() );
        addFieldData ( p, "p" );
    }

    // real and imaginary part
    scalars->SetNumberOfComponents ( 2 );

    // fill the scalars vector and add it to imageData_
    Teuchos::ArrayRCP<const std::complex<double> > xView;
    for ( int vec=0; vec<numVecs; vec++ )
    {
        xView = x.getVector ( vec )->get1dView();
        scalars->SetName ( scNames[vec].c_str() );
        for ( int k=0; k<numPoints; k++ )
        {
            if ( isScrambled )
            {
                // TODO replace by InsertNextTuple
                scalars->InsertNextValue ( p[k]>=0 ? std::real ( xView[p[k]] ) : dummy );
                scalars->InsertNextValue ( p[k]>=0 ? std::imag ( xView[p[k]] ) : dummy );
            }
            else
            {
                scalars->InsertNextValue ( std::real ( xView[k] ) );
                scalars->InsertNextValue ( std::imag ( xView[k] ) );
            }
        }
        imageData_->GetPointData()->AddArray ( scalars );
    }

    return;
}
// =============================================================================
