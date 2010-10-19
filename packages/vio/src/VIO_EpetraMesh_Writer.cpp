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

#include "VIO_EpetraMesh_Writer.h"

#include "VIO_EpetraMesh_Mesh.h"

#include <EpetraExt_Utils.h> // to_string
#include <Tpetra_Vector.hpp>
#include <Teuchos_Comm.hpp>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkTriangle.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

// =============================================================================
VIO::EpetraMesh::Writer::
Writer ( const std::string & filePath ) :
           VIO::Writer::Abstract( filePath ),
           mesh_( Teuchos::null ),
           comm_( Teuchos::null )
{
  vtkDataSet_ = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // analyze the file name for extension
  int         dotPos    = filePath.rfind ( "." );
  std::string extension = filePath.substr ( dotPos+1, filePath.size()-dotPos-1 );
  // convert to lower case
  std::transform( extension.begin(), extension.end(),
                  extension.begin(), ::tolower
                );

  if ( extension.compare ( "vtk" ) == 0 )
  {
  }
  else if ( extension.compare ( "vtu" ) == 0 )
  {
  }
  else if ( extension.compare ( "pvtu" ) == 0 )
  {
  }
  else
  {
      TEST_FOR_EXCEPTION ( true,
                            std::logic_error,
                            "Error when writing file \"" << filePath
                            << "\". File name extension \"" << extension << "\" "
                            << "not recognized. Must be one of \"vtk\", "
                            << "\"vtu\"." );
  }

  return;
}
// =============================================================================
VIO::EpetraMesh::Writer::
~Writer ()
{
}
// =============================================================================
void
VIO::EpetraMesh::Writer::
setMesh( const VIO::EpetraMesh::Mesh & mesh )
{
  mesh_ = Teuchos::rcp( new VIO::EpetraMesh::Mesh( mesh ) );

  // cast into a vtkUnstructuredGrid
  vtkSmartPointer<vtkUnstructuredGrid> vtkMesh =
      dynamic_cast<vtkUnstructuredGrid*> ( vtkDataSet_.GetPointer() );
  TEUCHOS_ASSERT_INEQUALITY( 0, !=, vtkMesh );

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
  vtkMesh->SetPoints( points );
  // ---------------------------------------------------------------------------
  // set cells
  vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();

//         MeshBase::const_element_iterator el     = mesh.elements_begin();
//   const MeshBase::const_element_iterator end_el = mesh.elements_end();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> > elems = mesh.getElems();
  Teuchos::ArrayRCP<const VIO::EpetraMesh::Mesh::ElementType> elemTypes = mesh.getElemTypes();
  for ( int k=0; k<elems.size(); k++ )
  {
      vtkIdType celltype = VTK_EMPTY_CELL; // initialize to something to avoid compiler warning

      switch( elemTypes[k] )
      {
              case VIO::EpetraMesh::Mesh::EDGE2:
                      celltype = VTK_LINE;
                      break;
              case VIO::EpetraMesh::Mesh::EDGE3:
                      celltype = VTK_QUADRATIC_EDGE;
                      break;// 1
              case VIO::EpetraMesh::Mesh::TRI3:
                      celltype = VTK_TRIANGLE;
                      break;// 3
              case VIO::EpetraMesh::Mesh::TRI6:
                      celltype = VTK_QUADRATIC_TRIANGLE;
                      break;// 4
              case VIO::EpetraMesh::Mesh::QUAD4:
                      celltype = VTK_QUAD;
                      break;// 5
              case VIO::EpetraMesh::Mesh::QUAD8:
                      celltype = VTK_QUADRATIC_QUAD;
                      break;// 6
              case VIO::EpetraMesh::Mesh::TET4:
                      celltype = VTK_TETRA;
                      break;// 8
              case VIO::EpetraMesh::Mesh::TET10:
                      celltype = VTK_QUADRATIC_TETRA;
                      break;// 9
              case VIO::EpetraMesh::Mesh::HEX8:
                      celltype = VTK_HEXAHEDRON;
                      break;// 10
              case VIO::EpetraMesh::Mesh::HEX20:
                      celltype = VTK_QUADRATIC_HEXAHEDRON;
                      break;// 12
              case VIO::EpetraMesh::Mesh::PRISM6:
                      celltype = VTK_WEDGE;
                      break;// 13
              case VIO::EpetraMesh::Mesh::PRISM15:
                      celltype = VTK_HIGHER_ORDER_WEDGE;
                      break;// 14
              case VIO::EpetraMesh::Mesh::PYRAMID5:
                      celltype = VTK_PYRAMID;
                      break;// 16
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
              case VIO::EpetraMesh::Mesh::QUAD9:
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

      vtkMesh->InsertNextCell( celltype, pts.GetPointer() );
  }
  // ---------------------------------------------------------------------------
  return;
}
// =============================================================================
void
VIO::EpetraMesh::Writer::
setValues( const Epetra_MultiVector          & x,
           const Teuchos::Array<std::string> & scalarsNames
         )
{
  unsigned int numVecs = x.NumVectors();
  unsigned int numPoints = x.GlobalLength();

  // cast into a vtkUnstructuredGrid
  vtkSmartPointer<vtkUnstructuredGrid> vtkMesh =
      dynamic_cast<vtkUnstructuredGrid*> ( vtkDataSet_.GetPointer() );
  TEUCHOS_ASSERT_INEQUALITY( 0, !=, vtkMesh );

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
      vtkMesh->GetPointData()->AddArray ( scalars );
  }

  return;
}
// =============================================================================
void
VIO::EpetraMesh::Writer::
setValues( const Tpetra::MultiVector<double> & x,
           const Teuchos::Array<std::string> & scalarsNames )
{
  unsigned int numVecs = x.getNumVectors();
  unsigned int numPoints = x.getGlobalLength();

  // cast into a vtkUnstructuredGrid
  vtkSmartPointer<vtkUnstructuredGrid> vtkMesh =
      dynamic_cast<vtkUnstructuredGrid*> ( vtkDataSet_.GetPointer() );
  TEUCHOS_ASSERT_INEQUALITY( 0, !=, vtkMesh );

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
      vtkMesh->GetPointData()->AddArray ( scalars );
  }

  return;
}
// =============================================================================
void
VIO::EpetraMesh::Writer::
setValues( const ComplexMultiVector          & z,
           const Teuchos::Array<std::string> & scalarsNames
         )
{
  unsigned int numVecs = z.getNumVectors();

  // cast into a vtkUnstructuredGrid
  vtkSmartPointer<vtkUnstructuredGrid> vtkMesh =
      dynamic_cast<vtkUnstructuredGrid*> ( vtkDataSet_.GetPointer() );
  TEUCHOS_ASSERT_INEQUALITY( 0, !=, vtkMesh );

  comm_ = z.getMap()->getComm();

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

  // real and imaginary part
  scalars->SetNumberOfComponents ( 2 );

  for ( int vec=0; vec<numVecs; vec++ )
  {
      Teuchos::ArrayRCP<const std::complex<double> > zView =
          z.getVector(vec)->get1dView();

      // fill the array
      for ( int k=0; k<zView.size(); k++ )
      {
          scalars->InsertNextValue ( zView[k].real() );
          scalars->InsertNextValue ( zView[k].imag() );
      }

      scalars->SetName ( scNames[vec].c_str() );

//       scalars->Print( std::cout );

      vtkMesh->GetPointData()->AddArray ( scalars );
  }

  return;
}
// =============================================================================
void
VIO::EpetraMesh::Writer::
write () const
{
    // cast mesh into a vtkUnstructuredGrid
    vtkSmartPointer<vtkUnstructuredGrid> vtkMesh =
        dynamic_cast<vtkUnstructuredGrid*> ( vtkDataSet_.GetPointer() );
    TEUCHOS_ASSERT_INEQUALITY( 0, !=, vtkMesh );


    TEUCHOS_ASSERT( !comm_.is_null() );

    // write the file
    if ( comm_->getSize() == 1 )
    {
        // serial writer
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer
                = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName ( filePath_.c_str() );
        writer->SetInput ( &*vtkMesh );
        writer->Write();
   }
   else
   {
       // parallel writer
       vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer
               = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
       //     if ( comm_->getRank() == 0 )
       //     {
       //       writer->SetNumberOfPieces( comm_->getSize() );
       //     }
       //     else
       //     {
       writer->SetNumberOfPieces( comm_->getSize() );
       //       writer->SetCompressor( 0 );
       //       writer->SetDataModeToAscii();
       //       writer->SetDebug( 'a' );
       writer->SetStartPiece( comm_->getRank() );
       writer->SetEndPiece( comm_->getRank() );
       //     }
       writer->SetFileName ( filePath_.c_str() );
       writer->SetInput ( &*vtkMesh );
       writer->Write();
   }

    return;
}
// =============================================================================
