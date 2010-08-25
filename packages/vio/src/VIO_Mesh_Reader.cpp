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
#include <vtkXMLPUnstructuredGridReader.h>

#include <vtkCell.h>
#include <vtkFeatureEdges.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPointData.h>
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
read ( Teuchos::RCP<ComplexMultiVector>              & z,
       Teuchos::RCP<Mesh>                            & mesh,
       Teuchos::ParameterList                        & fieldData,
       const Teuchos::RCP<const Teuchos::Comm<int> > & TComm
     )
{
  
  if ( TComm->getRank() == 1 ) std::cout << "a" << std::endl;

    // check the filename extension
    int         dotPos    = filePath_.rfind ( "." );
    std::string extension = filePath_.substr ( dotPos+1, filePath_.size()-dotPos-1 );
  
  if ( TComm->getRank() == 1 ) std::cout << "b" << std::endl;
  
    vtkSmartPointer<vtkUnstructuredGrid> vtkMesh;
  
  if ( TComm->getRank() == 1 ) std::cout << "c" << std::endl;
    
    if ( extension.compare ( "vtk" ) == 0 )
    {
        vtkSmartPointer<vtkUnstructuredGridReader> reader =
            vtkSmartPointer<vtkUnstructuredGridReader>::New();
        reader->SetFileName ( filePath_.c_str() );
        reader->Update();
        vtkMesh = reader->GetOutput();

    }
    else if ( extension.compare ( "vtu" ) == 0 )
    {
      
      if ( TComm->getRank() == 1 ) std::cout << "d" << std::endl;
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            if ( TComm->getRank() == 1 ) std::cout << "e" << std::endl;
        reader->SetFileName ( filePath_.c_str() );
      if ( TComm->getRank() == 1 ) std::cout << "f" << std::endl;
        reader->Update();
      if ( TComm->getRank() == 1 ) std::cout << "g" << std::endl;
        vtkMesh = reader->GetOutput();
      if ( TComm->getRank() == 1 ) std::cout << "h" << std::endl;
    }
    else if ( extension.compare ( "pvtu" ) == 0 )
    {
        vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
        reader->SetFileName ( filePath_.c_str() );
        reader->Update();
        vtkMesh = reader->GetOutput();
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
    
    
    if ( TComm->getRank() == 1 ) std::cout << "i" << std::endl;
    // read the data
    z         = this->extractStateData_ ( vtkMesh, TComm );
    if ( TComm->getRank() == 1 ) std::cout << "j" << std::endl;
    mesh      = this->extractMeshData_ ( vtkMesh, TComm );
    if ( TComm->getRank() == 1 ) std::cout << "k" << std::endl;
    fieldData = this->readFieldData_ ( vtkMesh );
    if ( TComm->getRank() == 1 ) std::cout << "l" << std::endl;
    
    return;
}
// =============================================================================
Teuchos::RCP<VIO::Mesh::Mesh>
VIO::Mesh::Reader::
extractMeshData_( const vtkSmartPointer<vtkUnstructuredGrid>    & vtkMesh,
                  const Teuchos::RCP<const Teuchos::Comm<int> > & TComm
                ) const
{ 
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp( new Mesh( TComm ) );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  // get points
  int numPoints = vtkMesh->GetNumberOfPoints();
  Teuchos::ArrayRCP<Point> points( numPoints );
  for ( unsigned int k=0; k<numPoints; k++ )
      vtkMesh->GetPoint( k, points[k].getRawPtr() );

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
  int globalNumElems = vtkMesh->GetNumberOfCells();

  // create an appropriate map
  Teuchos::RCP<const Tpetra::Map<ORD> > elemsMap =
      Teuchos::rcp ( new Tpetra::Map<ORD> ( globalNumElems, 0, TComm ) );
  
  int localNumElems = elemsMap->getNodeNumElements();
        
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<ORD> > elems( localNumElems );
  Teuchos::ArrayRCP<Mesh::ElementType>       elemTypes( localNumElems );
  
  for ( unsigned int k=0; k<localNumElems; k++ )
  {
      // set the connectivity table
      vtkCell * cell = vtkMesh->GetCell( elemsMap->getGlobalElement(k) );
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
  return mesh;
}
// =============================================================================
Teuchos::RCP<ComplexMultiVector>
VIO::Mesh::Reader::
extractStateData_ ( const vtkSmartPointer<vtkDataSet>             & vtkData,
                    const Teuchos::RCP<const Teuchos::Comm<int> > & TComm
                  ) const
{ 
    vtkIdType numArrays = vtkData->GetPointData()->GetNumberOfArrays();

    TEUCHOS_ASSERT_EQUALITY ( numArrays, 1 );

    const vtkSmartPointer<vtkDataArray> & array = 
        vtkData->GetPointData()->GetArray(0);
    
    vtkIdType numComponents = array->GetNumberOfComponents();

    TEUCHOS_ASSERT_EQUALITY ( numComponents, 2 );    // for *complex* values
    
    // this is the total number of grid points
    vtkIdType numPoints = array->GetNumberOfTuples();

    // create an appropriate map
    Teuchos::RCP<const Tpetra::Map<ORD> > ComplexMap =
        Teuchos::rcp ( new Tpetra::Map<ORD> ( numPoints, 0, TComm ) );

    Teuchos::RCP<ComplexMultiVector> z =
        Teuchos::rcp ( new ComplexVector ( ComplexMap ) );

    // fill z
    Teuchos::ArrayRCP<std::complex<double> > zLocalView = z->get1dViewNonConst();
    double val[2];
    for ( int k = 0; k < zLocalView.size(); k++ )
    {
        array->GetTuple( ComplexMap->getGlobalElement(k), val );
        zLocalView[k] = double_complex( val[0], val[1] );
    }

    return z;
}
// =============================================================================
// non-member function
void
VIO::Mesh::
read ( const Teuchos::RCP<const Teuchos::Comm<int> > & TComm,
       const std::string                             & filePath,
       Teuchos::RCP<ComplexMultiVector>              & z,
       Teuchos::RCP<VIO::Mesh::Mesh>                 & mesh,
       Teuchos::ParameterList                        & fieldData
     )
{
  VIO::Mesh::Reader reader( filePath );
  reader.read( z, mesh, fieldData, TComm );
  return;
}
// =============================================================================