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

#include "VIO_Mesh_Reader.h"

#include "VIO_Mesh_Mesh.h"

#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkCell.h>
#include <vtkFeatureEdges.h>
#include <vtkDataSetSurfaceFilter.h>

#include <vtkCellArray.h>

// =============================================================================
VIO::Mesh::Reader::
Reader ( const std::string & filePath ):
  VIO::Reader::Abstract( filePath )
{
}
// =============================================================================
VIO::Mesh::Reader::
~Reader ()
{
}
// =============================================================================
void
VIO::Mesh::Reader::
read ( Teuchos::RCP<Mesh> & mesh
     )
{
    // check the filename extension
    int         dotPos    = filePath_.rfind ( "." );
    std::string extension = filePath_.substr ( dotPos+1, filePath_.size()-dotPos-1 );
  
    vtkSmartPointer<vtkUnstructuredGrid> vtkGrid;
    
    if ( extension.compare ( "vtk" ) == 0 )
    {
        vtkSmartPointer<vtkUnstructuredGridReader> reader =
            vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName ( filePath_.c_str() );
        reader->Update();
        vtkGrid = reader->GetOutput();

    }
    else if ( extension.compare ( "vtu" ) == 0 )
    {
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader->SetFileName ( filePath_.c_str() );
        reader->Update();
        vtkGrid = reader->GetOutput();
    }
    else
    {
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Error when reading file \"" << filePath_
                             << "\". File name extension \"" << extension << "\" "
                             << "not recognized. Must be one of \"vtk\", "
                             << "\"vtu\"." );
    }
    
    this->processVtkData_ ( vtkGrid,
                            mesh );

    return;
}
// =============================================================================
void
VIO::Mesh::Reader::
processVtkData_( const vtkSmartPointer<vtkUnstructuredGrid> & vtkMesh,
                 Teuchos::RCP<Mesh>                         & mesh
               )
{ 
  TEUCHOS_ASSERT( !mesh.is_null() );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  // get points
  int numPoints = vtkMesh->GetNumberOfPoints();
  Teuchos::ArrayRCP<Point> points( numPoints );
  for ( unsigned int k=0; k<numPoints; k++ )
  {
      vtkMesh->GetPoint( k, points[k].getRawPtr() );
  }
  mesh->setNodes( points );
  
//   for ( int k=0; k<points.size(); k++ )
//       std::cout << points[k] << std::endl;
//   
//   std::cout << "VVV" << std::endl;
//   vtkMesh->Print( std::cout );
  
  // determine the boundary points
  // transform vtkMesh into vtkPolyData
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = 
      vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfaceFilter->SetInput(vtkMesh);
  surfaceFilter->Update(); 
//   vtkPolyData* polydata = surfaceFilter->GetOutput();

  // filter out the boundary edges
  vtkFeatureEdges * pEdges = vtkFeatureEdges::New();
  pEdges->SetInput( surfaceFilter->GetOutput() );
  pEdges->BoundaryEdgesOn();
  pEdges->FeatureEdgesOff();
  pEdges->NonManifoldEdgesOff();
  pEdges->ManifoldEdgesOff();
  pEdges->Update();
  
  vtkPolyData * poly = pEdges->GetOutput();
  
//   pEdges->ColoringOff();
// //   pEdges->GetOutput()->BuildCells();
//   
//   std::cout << "V1" << std::endl;
//   std::cout << pEdges->GetOutput()->GetNumberOfCells() << std::endl;
//   std::cout << pEdges->GetOutput()->GetNumberOfLines() << std::endl;
//   std::cout << pEdges->GetOutput()->GetNumberOfPieces() << std::endl;
//   std::cout << pEdges->GetOutput()->GetNumberOfPolys() << std::endl;
//   std::cout << pEdges->GetOutput()->GetNumberOfStrips() << std::endl;
//   std::cout << pEdges->GetOutput()->GetNumberOfVerts() << std::endl;
//   std::cout << pEdges->GetOutput()->GetNumberOfPoints() << std::endl;
// //   std::cout << pEdges->GetOutput()->GetNumberOfElements() << std::endl;
//   
// //   pEdges->GetOutput()->GetLines()->GetCell(0)->Print( std::cout );
//   
// //   pEdges->GetOutput()->BuildCells();
//   TEUCHOS_ASSERT_EQUALITY( 0, pEdges->GetErrorCode() );
// 
//   std::cout << "W-1" << std::endl;
//   int numLines = poly->GetNumberOfLines();
//   vtkIdType * pts;
// //   vtkIdList * pts2;
//   for ( vtkIdType lineId=0; lineId<numLines; lineId++ )
//   {
//       vtkIdType numPoints;
//       pEdges->GetOutput()->GetLines()->GetCell( lineId, numPoints, pts );
// //       pEdges->GetOutput()->GetCellPoints( cellId, pts2 );
//       
//       if ( numPoints!=2 )
//       {
//         std::cout << "AAAAAAAAH!" << std::endl;
//         continue;
//       }
//       // make sure we're dealing with an actual *edge* here.
//       TEUCHOS_ASSERT_EQUALITY( numPoints, 2 );
//       std::cout << pts[0] << " " << pts[1] << std::endl;
//   }
  
//   std::cout << "W0" << std::endl;
//   poly->GetLines()->Print( std::cout );
//   
//   std::cout << "W1" << std::endl;
//   poly->GetPoints()->Print( std::cout );
//   
//   std::cout << "W2" << std::endl;

  double x[3];
//   Teuchos::ArrayRCP<bool> boundaryPoints( numPoints, false );
  Teuchos::ArrayRCP<bool> isBoundaryNodes( numPoints );
  for ( int k=0; k<poly->GetNumberOfPoints(); k++ )
  {
      poly->GetPoint( k, x );
      vtkIdType ptId = vtkMesh->FindPoint( x );
      isBoundaryNodes[ ptId ] = true;
  }
  mesh->setBoundaryNodes( isBoundaryNodes );
  
//   for ( int k=0; k<boundaryPoints.size(); k++ )
//     std::cout << boundaryPoints[k] << std::endl;
  
//   std::cout << "WWW" << std::endl;
//   poly->Print( std::cout );
//   vtkCellArray * verts = poly->GetVerts();
//   verts->Print( std::cout );
  
//   std::cout << "XXX" << std::endl;
//   poly->GetData()->Print( std::cout );
  
//   int numBoundaryPoints = poly->GetNumberOfPoints();
//   Teuchos::ArrayRCP<Teuchos::Tuple<double,3> > boundaryPoints( numBoundaryPoints );
//   for ( unsigned int k=0; k<numBoundaryPoints; k++ )
//       poly->GetPoint( k, boundaryPoints[k].getRawPtr() );
//   mesh->setBoundaryNodes( boundaryPoints );
  
//   for ( int k=0; k<points.size(); k++ )
//       std::cout << boundaryPoints[k] << std::endl;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // get cells
  int numElems = vtkMesh->GetNumberOfCells();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems( numElems );
  Teuchos::ArrayRCP<Mesh::ElementType>       elemTypes( numElems );
  for ( unsigned int k=0; k<numElems; k++ )
  {
      // set the connectivity table
      vtkCell * cell = vtkMesh->GetCell( k );
      int numPoints = cell->GetNumberOfPoints();
      elems[k] = Teuchos::ArrayRCP<ORD>( numPoints );
      for ( unsigned int l=0; l<numPoints; l++ )
          elems[k][l] = cell->GetPointId( l );

      // set the element type
      switch( cell->GetCellType() )
      {
          case VTK_LINE:
                  elemTypes[k] = Mesh::EDGE2;
                  break;
          case VTK_QUADRATIC_EDGE:
                  elemTypes[k] = Mesh::EDGE3;
                  break;
          case VTK_TRIANGLE:
                  elemTypes[k] = Mesh::TRI3;
                  break;
          case VTK_QUADRATIC_TRIANGLE:
                  elemTypes[k] = Mesh::TRI6;
                  break;
          case VTK_QUAD:
                  elemTypes[k] = Mesh::QUAD4;
                  break;
          case VTK_QUADRATIC_QUAD:
                  elemTypes[k] = Mesh::QUAD8;
                  break;
          case VTK_BIQUADRATIC_QUAD:
                  elemTypes[k] = Mesh::QUAD9;
                  break;
          case VTK_TETRA:
                  elemTypes[k] = Mesh::TET4;
                  break;
          case VTK_QUADRATIC_TETRA:
                  elemTypes[k] = Mesh::TET10;
                  break;
          case VTK_HEXAHEDRON:
                  elemTypes[k] = Mesh::HEX8;
                  break;
          case VTK_QUADRATIC_HEXAHEDRON:
                  elemTypes[k] = Mesh::HEX20;
                  break;
          case VTK_WEDGE:
                  elemTypes[k] = Mesh::PRISM6;
                  break;
          case VTK_HIGHER_ORDER_WEDGE:
                  elemTypes[k] = Mesh::PRISM15;
                  break;
          case VTK_PYRAMID:
                  elemTypes[k] = Mesh::PYRAMID5;
                  break;
          default:
              TEST_FOR_EXCEPTION( true,
                                  std::logic_error,
                                  "Unknown type \""<< cell->GetCellType()  <<"\"." );
      }
  }
  mesh->setElems( elems );
  mesh->setElemTypes( elemTypes );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return;
}
// =============================================================================