/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

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

#include "VIO_Mesh_Writer.h"

#include "VIO_Mesh_Mesh.h"

#include <EpetraExt_Utils.h> // to_string
#include <Tpetra_Vector.hpp>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkTriangle.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

// =============================================================================
VIO::Mesh::Writer::
Writer ( const std::string & filePath ) :
           filePath_ ( filePath ),
           vtkMesh_( vtkSmartPointer<vtkUnstructuredGrid>::New() ),
           mesh_( Teuchos::null )
{
}
// =============================================================================
VIO::Mesh::Writer::
~Writer ()
{
}
// =============================================================================
void
VIO::Mesh::Writer::
setMesh( const VIO::Mesh::Mesh & mesh )
{
  mesh_ = Teuchos::rcp( new VIO::Mesh::Mesh( mesh ) );
  
  // ---------------------------------------------------------------------------
  // set points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//   std::cout << "\n\nPoint order in VTK:" << std::endl;
  Teuchos::ArrayRCP<Point> nodes = mesh.getNodes();
  for ( int k=0; k<nodes.size(); k++ )
  {
      points->InsertNextPoint( nodes[k].getRawPtr()  );
//       std::cout << nodes[k] << std::endl;
  }
  vtkMesh_->SetPoints( points );
  // ---------------------------------------------------------------------------
  // set cells
  vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();

//         MeshBase::const_element_iterator el     = mesh.elements_begin();
//   const MeshBase::const_element_iterator end_el = mesh.elements_end();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> >   elems = mesh.getElems();
  Teuchos::ArrayRCP<const VIO::Mesh::Mesh::ElementType> elemTypes = mesh.getElemTypes();
  for ( int k=0; k<elems.size(); k++ )
  {
      vtkIdType celltype = VTK_EMPTY_CELL; // initialize to something to avoid compiler warning

      switch( elemTypes[k] )
      {
              case VIO::Mesh::Mesh::EDGE2:
                      celltype = VTK_LINE;
                      break;
              case VIO::Mesh::Mesh::EDGE3:
                      celltype = VTK_QUADRATIC_EDGE;
                      break;// 1
              case VIO::Mesh::Mesh::TRI3:
                      celltype = VTK_TRIANGLE;
                      break;// 3
              case VIO::Mesh::Mesh::TRI6:
                      celltype = VTK_QUADRATIC_TRIANGLE;
                      break;// 4
              case VIO::Mesh::Mesh::QUAD4:
                      celltype = VTK_QUAD;
                      break;// 5
              case VIO::Mesh::Mesh::QUAD8:
                      celltype = VTK_QUADRATIC_QUAD;
                      break;// 6
              case VIO::Mesh::Mesh::TET4:
                      celltype = VTK_TETRA;
                      break;// 8
              case VIO::Mesh::Mesh::TET10:
                      celltype = VTK_QUADRATIC_TETRA;
                      break;// 9
              case VIO::Mesh::Mesh::HEX8:
                      celltype = VTK_HEXAHEDRON;
                      break;// 10
              case VIO::Mesh::Mesh::HEX20:
                      celltype = VTK_QUADRATIC_HEXAHEDRON;
                      break;// 12
              case VIO::Mesh::Mesh::PRISM6:
                      celltype = VTK_WEDGE;
                      break;// 13
              case VIO::Mesh::Mesh::PRISM15:
                      celltype = VTK_HIGHER_ORDER_WEDGE;
                      break;// 14
              case VIO::Mesh::Mesh::PYRAMID5:
                      celltype = VTK_PYRAMID;
                      break;// 16
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
              case VIO::Mesh::Mesh::QUAD9:
                      celltype = VTK_BIQUADRATIC_QUAD;
                      break;
#else
              case QUAD9:
#endif
              default:
                  TEST_FOR_EXCEPTION( true,
                                      std::logic_error,
                                      "Element type \""<< elemTypes[k] <<"\" not implemented." );
      }

      pts->SetNumberOfIds( elems[k].size() );

      // get the connectivity for this element
      for(unsigned int i=0; i<elems[k].size(); i++ )
          pts->InsertId( i, elems[k][i] );

      vtkMesh_->InsertNextCell( celltype, pts.GetPointer() );
  }
  // ---------------------------------------------------------------------------
  return; 
}
// =============================================================================
void
VIO::Mesh::Writer::
setValues( const Epetra_MultiVector          & x,
           const Teuchos::Array<std::string> & scalarsNames )
{
  unsigned int numVecs = x.NumVectors();
  unsigned int numPoints = x.GlobalLength();
  
  // get scalarsNames, and insert default names if empty
  Teuchos::Array<std::string> scNames ( scalarsNames );
  if ( scNames.empty() )
  {
      scNames.resize ( numVecs );
      for ( int vec=0; vec<numVecs; vec++ )
          scNames[vec] = "x" + EpetraExt::toString ( vec );
  }
  
  // fill the scalar field
  vtkSmartPointer<vtkDoubleArray> scalars =
      vtkSmartPointer<vtkDoubleArray>::New();

  for ( int vec=0; vec<numVecs; vec++ )
  { 
      scalars->SetName ( scNames[vec].c_str() );
      for ( int k=0; k<numPoints; k++ )
      {
//           const unsigned int dof_id = libmeshMesh_->node(k).dof_number(0,k,0);
          scalars->InsertNextValue ( x[vec][k] );
      }
      vtkMesh_->GetPointData()->AddArray ( scalars );
  }
  
  return;
}
// =============================================================================
void
VIO::Mesh::Writer::
setValues( const Tpetra::MultiVector<double> & x,
           const Teuchos::Array<std::string> & scalarsNames )
{
  unsigned int numVecs = x.getNumVectors();
  unsigned int numPoints = x.getGlobalLength();
  
  // get scalarsNames, and insert default names if empty
  Teuchos::Array<std::string> scNames ( scalarsNames );
  if ( scNames.empty() )
  {
      scNames.resize ( numVecs );
      for ( int vec=0; vec<numVecs; vec++ )
          scNames[vec] = "x" + EpetraExt::toString ( vec );
  }
  
  // fill the scalar field
  vtkSmartPointer<vtkDoubleArray> scalars =
      vtkSmartPointer<vtkDoubleArray>::New();

  for ( int vec=0; vec<numVecs; vec++ )
  {
      Teuchos::ArrayRCP<const double> xView = x.getVector(vec)->get1dView(); 
      scalars->SetName ( scNames[vec].c_str() );
      for ( int k=0; k<numPoints; k++ )
      {
//           const unsigned int dof_id = libmeshMesh_->node(k).dof_number(0,k,0);
          scalars->InsertNextValue ( xView[k] );
      }
      vtkMesh_->GetPointData()->AddArray ( scalars );
  }
  
  return;
}
// =============================================================================
void
VIO::Mesh::Writer::
setValues( const ComplexMultiVector          & z,
           const Teuchos::Array<std::string> & scalarsNames
         )
{
  unsigned int numVecs = z.getNumVectors();
  unsigned int numPoints = z.getGlobalLength();
  
  // get scalarsNames, and insert default names if empty
  Teuchos::Array<std::string> scNames ( scalarsNames );
  if ( scNames.empty() )
  {
      scNames.resize ( numVecs );
      for ( int vec=0; vec<numVecs; vec++ )
          scNames[vec] = "z" + EpetraExt::toString ( vec );
  }
  
  // fill the scalar field
  vtkSmartPointer<vtkDoubleArray> scalars =
      vtkSmartPointer<vtkDoubleArray>::New();

  std::stringstream strstream;
  for ( int vec=0; vec<numVecs; vec++ )
  {
      Teuchos::ArrayRCP<const std::complex<double> > zView =
          z.getVector(vec)->get1dView();
          
      // real part
      strstream << scNames[vec] << "_real";
      scalars->SetName ( strstream.str().c_str() );
      strstream.flush();
      for ( int k=0; k<numPoints; k++ )
          scalars->InsertNextValue ( zView[k].real() );
      vtkMesh_->GetPointData()->AddArray ( scalars );
      
      // imaginary part
      strstream << scNames[vec] << "_imag";
      scalars->SetName ( strstream.str().c_str() );
      strstream.flush();
      for ( int k=0; k<numPoints; k++ )
          scalars->InsertNextValue ( zView[k].imag() );
      vtkMesh_->GetPointData()->AddArray ( scalars );
  }
  
  return;
}
// =============================================================================
// void
// VIO::Mesh::Writer::
// addFieldData ( const Teuchos::Array<int> & array,
//                const std::string         & name )
// {
//     // create field data
//     vtkSmartPointer<vtkIntArray> fieldData = vtkSmartPointer<vtkIntArray>::New();
// 
//     fieldData->SetName ( name.c_str() );
// 
//     // fill the field
//     for ( int k=0; k<array.length(); k++ )
//         fieldData->InsertNextValue ( array[k] );
// 
//     imageData_->GetFieldData()->AddArray ( fieldData );
// 
//     return;
// }
// // =============================================================================
// void
// VIO::Mesh::Writer::
// addParameterList ( const Teuchos::ParameterList & problemParams )
// {
//     // add to imageData_
//     Teuchos::map<std::string, Teuchos::ParameterEntry>::const_iterator k;
//     for ( k = problemParams.begin(); k != problemParams.end(); ++k )
//     {
//         std::string paramName = problemParams.name ( k );
//         if ( problemParams.isType<int> ( paramName ) )
//         {
//             vtkSmartPointer<vtkIntArray> fieldData = vtkSmartPointer<vtkIntArray>::New();
//             fieldData->InsertNextValue ( problemParams.get<int> ( paramName ) );
//             fieldData->SetName ( paramName.c_str() );
//             imageData_->GetFieldData()->AddArray ( fieldData );
//         }
//         else if ( problemParams.isType<double> ( paramName ) )
//         {
//             vtkSmartPointer<vtkDoubleArray> fieldData = vtkSmartPointer<vtkDoubleArray>::New();
//             fieldData->InsertNextValue ( problemParams.get<double> ( paramName ) );
//             fieldData->SetName ( paramName.c_str() );
//             imageData_->GetFieldData()->AddArray ( fieldData );
//         }
//         else if ( problemParams.isType<Teuchos::Array<double> > ( paramName ) )
//         {
//             vtkSmartPointer<vtkDoubleArray> fieldData = vtkSmartPointer<vtkDoubleArray>::New();
//             const Teuchos::Array<double> & arr = problemParams.get<Teuchos::Array<double> > ( paramName );
//             for ( int k=0; k<arr.length(); k++ )
//                 fieldData->InsertNextValue ( arr[k] );
//             fieldData->SetName ( paramName.c_str() );
//             imageData_->GetFieldData()->AddArray ( fieldData );
//         }
//         else
//         {
//             TEST_FOR_EXCEPTION ( true,
//                                  std::runtime_error,
//                                  "Illegal type of parameter \"" << paramName << "\"." );
//         }
//     }
// 
//     return;
// }
// // =============================================================================
// void
// VIO::Mesh::Writer::
// setImageData ( const Epetra_MultiVector              & x,
//                const Teuchos::Tuple<unsigned int,2>  & Nx,
//                const Teuchos::Tuple<double,2>        & h,
//                const Teuchos::Array<int>             & p,
//                const Teuchos::Array<std::string>     & scalarsNames
//              )
// {
//     int numVecs   = x.NumVectors();
//     int numPoints = ( Nx[0]+1 ) * ( Nx[1]+1 );
// 
//     // get scalarsNames, and insert default names if empty
//     Teuchos::Array<std::string> scNames ( scalarsNames );
//     if ( scNames.empty() )
//     {
//         scNames.resize ( numVecs );
//         for ( int vec=0; vec<numVecs; vec++ )
//             scNames[vec] = "x" + EpetraExt::toString ( vec );
//     }
// 
//     // set other image data
//     imageData_->SetDimensions ( Nx[0]+1, Nx[1]+1, 1 );
//     imageData_->SetOrigin ( 0.0, 0.0, 0.0 );
//     imageData_->SetSpacing ( h[0], h[1], 0.0 );
// 
//     // fill the scalar field
//     vtkSmartPointer<vtkDoubleArray> scalars =
//         vtkSmartPointer<vtkDoubleArray>::New();
// 
//     bool isScrambled = !p.empty();
// 
//     if ( isScrambled )
//     {
//         TEUCHOS_ASSERT_EQUALITY ( numPoints, p.length() );
//         addFieldData ( p, "p" );
//     }
// 
//     // fill the scalars vector and add it to imageData_
//     if ( isScrambled )
//     {
//         double dummy = 0.0;
//         for ( int vec=0; vec<numVecs; vec++ )
//         {
//             scalars->SetName ( scNames[vec].c_str() );
//             for ( int k=0; k<numPoints; k++ )
//                 scalars->InsertNextValue ( p[k]>=0 ? x[vec][p[k]] : dummy );
//             imageData_->GetPointData()->AddArray ( scalars );
//         }
//     }
//     else
//         for ( int vec=0; vec<numVecs; vec++ )
//         {
//             scalars->SetName ( scNames[vec].c_str() );
//             for ( int k=0; k<numPoints; k++ )
//                 scalars->InsertNextValue ( x[vec][k] );
//             imageData_->GetPointData()->AddArray ( scalars );
//         }
// 
//     return;
// }
// // =============================================================================
// void
// VIO::Mesh::Writer::
// setImageData ( const DoubleMultiVector               & x,
//                const Teuchos::Tuple<unsigned int,2>  & Nx,
//                const Teuchos::Tuple<double,2>        & h,
//                const Teuchos::Array<int>             & p,
//                const Teuchos::Array<std::string>     & scalarsNames
//              )
// {
//     int numVecs   = x.getNumVectors();
//     int numPoints = ( Nx[0]+1 ) * ( Nx[1]+1 );
// 
//     // get scalarsNames, and insert default names if empty
//     Teuchos::Array<std::string> scNames ( scalarsNames );
//     if ( scNames.empty() )
//     {
//         scNames.resize ( numVecs );
//         for ( int vec=0; vec<numVecs; vec++ )
//             scNames[vec] = "x" + EpetraExt::toString ( vec );
//     }
// 
//     // set other image data
//     imageData_->SetDimensions ( Nx[0]+1, Nx[1]+1, 1 );
//     imageData_->SetOrigin ( 0.0, 0.0, 0.0 );
//     imageData_->SetSpacing ( h[0], h[1], 0.0 );
// 
//     // fill the scalar field
//     vtkSmartPointer<vtkDoubleArray> scalars =
//         vtkSmartPointer<vtkDoubleArray>::New();
// 
//     double dummy = 0.0;
//     bool isScrambled = !p.empty();
//     if ( isScrambled )
//     {
//         TEUCHOS_ASSERT_EQUALITY ( numPoints, p.length() );
//         addFieldData ( p, "p" );
//     }
// 
//     // fill the scalars vector and add it to imageData_
//     Teuchos::ArrayRCP<const double> xView;
//     for ( int vec=0; vec<numVecs; vec++ )
//     {
//         xView = x.getVector ( vec )->get1dView();
//         scalars->SetName ( scNames[vec].c_str() );
//         for ( int k=0; k<numPoints; k++ )
//         {
//             if ( isScrambled )
//                 scalars->InsertNextValue ( p[k]>=0 ? xView[p[k]] : dummy );
//             else
//                 scalars->InsertNextValue ( xView[k] );
//         }
//         imageData_->GetPointData()->AddArray ( scalars );
//     }
// 
//     return;
// }
// // =============================================================================
// void
// VIO::Mesh::Writer::
// setImageData ( const ComplexMultiVector              & x,
//                const Teuchos::Tuple<unsigned int,2>  & Nx,
//                const Teuchos::Tuple<double,2>        & h,
//                const Teuchos::Array<int>             & p,
//                const Teuchos::Array<std::string>     & scalarsNames
//              )
// {
//     int numVecs   = x.getNumVectors();
//     int numPoints = ( Nx[0]+1 ) * ( Nx[1]+1 );
// 
//     // get scalarsNames, and insert default names if empty
//     Teuchos::Array<std::string> scNames ( scalarsNames );
//     if ( scNames.empty() )
//     {
//         scNames.resize ( numVecs );
//         for ( int vec=0; vec<numVecs; vec++ )
//             scNames[vec] = "z" + EpetraExt::toString ( vec );
//     }
// 
//     // set other image data
//     imageData_->SetDimensions ( Nx[0]+1, Nx[1]+1, 1 );
//     imageData_->SetOrigin ( 0.0, 0.0, 0.0 );
//     imageData_->SetSpacing ( h[0], h[1], 0.0 );
// 
//     // fill the scalar field
//     vtkSmartPointer<vtkDoubleArray> scalars =
//         vtkSmartPointer<vtkDoubleArray>::New();
// 
//     double dummy = 0.0;
//     bool isScrambled = !p.empty();
//     if ( isScrambled )
//     {
//         TEUCHOS_ASSERT_EQUALITY ( numPoints, p.length() );
//         addFieldData ( p, "p" );
//     }
// 
//     // real and imaginary part
//     scalars->SetNumberOfComponents ( 2 );
// 
//     // fill the scalars vector and add it to imageData_
//     Teuchos::ArrayRCP<const std::complex<double> > xView;
//     for ( int vec=0; vec<numVecs; vec++ )
//     {
//         xView = x.getVector ( vec )->get1dView();
//         scalars->SetName ( scNames[vec].c_str() );
//         for ( int k=0; k<numPoints; k++ )
//         {
//             if ( isScrambled )
//             {
//                 // TODO replace by InsertNextTuple
//                 scalars->InsertNextValue ( p[k]>=0 ? std::real ( xView[p[k]] ) : dummy );
//                 scalars->InsertNextValue ( p[k]>=0 ? std::imag ( xView[p[k]] ) : dummy );
//             }
//             else
//             {
//                 scalars->InsertNextValue ( std::real ( xView[k] ) );
//                 scalars->InsertNextValue ( std::imag ( xView[k] ) );
//             }
//         }
//         imageData_->GetPointData()->AddArray ( scalars );
//     }
// 
//     return;
// }
// =============================================================================
void
VIO::Mesh::Writer::
write () const
{
    // write the file
//     vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
//         vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    writer->SetFileName ( filePath_.c_str() );
    writer->SetInput ( &*vtkMesh_ );
    writer->Write();

    return;
}
// =============================================================================
