#include "stateFileWriter.h"

#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLParameterListWriter.hpp>
#include <EpetraExt_Utils.h> // for toString

// =============================================================================
// Constructor
StateFileWriter::StateFileWriter( StaggeredGrid &sGrid ):
  SGrid(sGrid) // make this &sGrid
{
}
// =============================================================================



// =============================================================================
// Destructor
StateFileWriter::~StateFileWriter()
{
}
// =============================================================================



// =============================================================================
void StateFileWriter::writeFile( const std::string                 &fileName,
                                 const std::string                 &fileFormat,
                                 const std::vector<double_complex> &psi,
                                 const Teuchos::ParameterList      &problemParams )
{
  if ( !fileFormat.compare("legacyVTK") )
      writeLegacyVtkFile( fileName, psi, problemParams );
  else if ( !fileFormat.compare("VTI") )
      writeVtiFile( fileName, psi, problemParams );
  else if ( !fileFormat.compare("XDMF") )
      writeXdmfFile( fileName, psi, problemParams );
  else {
      std::cerr << "StateFileWriter::writeFile :" << std::endl
                << "Illegal file format \"" << fileFormat <<"\"." << std::endl
                << "Use one of \"legacyVTK\", \"VTI\", and \"XDMF\"." << std::endl;
      exit(EXIT_FAILURE);
  }
}
// =============================================================================



// =============================================================================
void StateFileWriter::writeLegacyVtkFile( const std::string                 &filename,
                                          const std::vector<double_complex> &psi,
                                          const Teuchos::ParameterList      &problemParams )
{
  int           Nx = SGrid.getNx(),
                k,
                index[2];
  double        h  = SGrid.getH();
  std::ofstream vtkfile;

  // open the file
  vtkfile.open( filename.c_str() );

  // write the VTK header
  vtkfile << "# vtk DataFile Version 2.0\n"
          << "dummy line\n"
          << "ASCII\n"
          << "DATASET STRUCTURED_POINTS\n"
          << "DIMENSIONS " << Nx+1 << " " << Nx+1 << " " << 1 << "\n"
          << "ORIGIN 0 0 0\n"
          << "SPACING " << h << " " << h << " " << 0.0 << "\n"
          << "POINT_DATA " << SGrid.getNumComplexUnknowns() << "\n";

  // Note that, when writing the data, the values of psi are assumed to be
  // given in lexicographic ordering.

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write abs(psi)
  vtkfile << "SCALARS abs(psi) float\n"
          << "LOOKUP_TABLE default\n";
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = SGrid.i2k( index );
          vtkfile << abs(psi[k]) << "\n";
      }
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write arg(psi)
  vtkfile << "SCALARS arg(psi) float\n"
          << "LOOKUP_TABLE default\n";
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = SGrid.i2k( index );
          vtkfile << arg(psi[k]) << "\n";
      }
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // close the file
  vtkfile.close();
}
// =============================================================================




// =============================================================================
void StateFileWriter::writeVtiFile( const std::string                 &filename,
                                    const std::vector<double_complex> &psi,
                                    const Teuchos::ParameterList      &problemParams )
{
  int    Nx = SGrid.getNx(),
         k,
         index[2];
  double h  = SGrid.getH();

  std::string str;

  Teuchos::XMLObject xmlPointData("PointData");
  xmlPointData.addAttribute( "Scalars", "abs(psi)" );

  // first build the XML structure
  Teuchos::XMLObject xmlDataArrayAbs("DataArray");
  xmlDataArrayAbs.addAttribute( "type", "Float32" );
  xmlDataArrayAbs.addAttribute( "Name", "abs(psi)" );
  xmlDataArrayAbs.addAttribute( "format", "ascii" );
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = SGrid.i2k( index );
          xmlDataArrayAbs.addContent( EpetraExt::toString(abs(psi[k])) + " ");
      }
  }
  xmlPointData.addChild(xmlDataArrayAbs);

  Teuchos::XMLObject xmlDataArrayArg("DataArray");
  xmlDataArrayArg.addAttribute( "type", "Float32" );
  xmlDataArrayArg.addAttribute( "Name", "arg(psi)" );
  xmlDataArrayArg.addAttribute( "format", "ascii" );
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = SGrid.i2k( index );
          xmlDataArrayArg.addContent( EpetraExt::toString(arg(psi[k])) + " " );
      }
  }
  xmlPointData.addChild(xmlDataArrayArg);

  Teuchos::XMLObject xmlPiece("Piece");
  str = "0 " + EpetraExt::toString(Nx) + " 0 " + EpetraExt::toString(Nx) + " 0 0";
  xmlPiece.addAttribute( "Extent", str );
  xmlPiece.addChild(xmlPointData);

  Teuchos::XMLObject xmlImageData("ImageData");
  xmlImageData.addAttribute( "WholeExtent", str );
  xmlImageData.addAttribute( "Origin", "0 0 0" );
  str = EpetraExt::toString(h) + " " + EpetraExt::toString(h) + " 0";
  xmlImageData.addAttribute( "Spacing", str );
  xmlImageData.addChild(xmlPiece);


  // append the problem parameters in XML form to the file
  Teuchos::XMLObject xmlParameterList("");
  xmlParameterList = Teuchos::XMLParameterListWriter::XMLParameterListWriter()
                                                        .toXML( problemParams );

  // define top level object
  Teuchos::XMLObject vtuxml("VTKFile");

  // append the parameter list to the embracing VTK XML object
  vtuxml.addChild(xmlParameterList);

  vtuxml.addAttribute( "type", "ImageData" );
  vtuxml.addAttribute( "version", "0.1" );
  vtuxml.addAttribute( "byte_order", "LittleEndian" );
  vtuxml.addChild(xmlImageData);

  // ---------------------------------------------------------------------------
  // write the contents to the file
  // open the file
  std::ofstream  vtkfile;
  vtkfile.open( filename.c_str() );

  // Do not plot the XML header as Teuchos' XML reader can't deal with it
  // vtkfile << "<?xml version=\"1.0\"?>" << std::endl;

  vtkfile << vtuxml;
  // close the file
  vtkfile.close();
  // ---------------------------------------------------------------------------
}
// =============================================================================




// =============================================================================
void StateFileWriter::writeXdmfFile( const std::string                 &filename,
                                     const std::vector<double_complex> &psi,
                                     const Teuchos::ParameterList      &problemParams )
{

  std::cerr << "StateFileWriter::writeXdmfFile not yet implemented." << std::endl;
  exit(EXIT_FAILURE);


// TODO: use toString from EpetraExt instead of those nasty ofstreams!

  int    Nx = SGrid.getNx(),
         k,
         index[2];
  double h  = SGrid.getH();
  std::ostringstream os;
  std::ofstream      xdmfFile;

  // ---------------------------------------------------------------------------
  // write the XDMF file
  // ---------------------------------------------------------------------------
  // set grid topology
  Teuchos::XMLObject xmlTopology("Topology");
  xmlTopology.addAttribute( "TopologyType", "3DCORECTMESH" );
  os << "1 " << Nx+1 << " " << Nx+1;
  xmlTopology.addAttribute( "Dimensions", os.str() );
  os.str("");

  // define origin
  Teuchos::XMLObject xmlOrigin("DataItem");
  xmlOrigin.addAttribute( "Name", "Origin" );
  xmlOrigin.addAttribute( "NumberType", "Float" );
  xmlOrigin.addAttribute( "Dimensions", "3" );
  xmlOrigin.addAttribute( "Format", "XML" );
  xmlOrigin.addContent( "0 0 0" );

  // define spacing
  Teuchos::XMLObject xmlSpacing("DataItem");
  xmlSpacing.addAttribute( "Name", "Spacing" );
  xmlSpacing.addAttribute( "NumberType", "Float" );
  xmlSpacing.addAttribute( "Dimensions", "3" );
  xmlSpacing.addAttribute( "Format", "XML" );
  os << "0 " << h << " " << h;
  xmlSpacing.addContent( os.str() );
  os.str("");

  // merge the latter two into geometry
  Teuchos::XMLObject xmlGeometry("Geometry");
  xmlGeometry.addAttribute( "Type", "ORIGIN_DXDYDZ" );
  xmlGeometry.addChild( xmlOrigin );
  xmlGeometry.addChild( xmlSpacing );

  // tell me where the actual ABS(PSI) data sits
  Teuchos::XMLObject xmlAbsData("DataItem");
  xmlAbsData.addAttribute( "NumberType", "Float" );
  xmlAbsData.addAttribute( "Precision", "4" );
  os << "1 " << Nx+1 << " " << Nx+1;
  xmlAbsData.addAttribute( "Dimensions", os.str() );
  os.str("");
  xmlAbsData.addAttribute( "Format", "HDF" );
  xmlAbsData.addContent( "myfile.h5:/abs(psi)" );

  Teuchos::XMLObject xmlAbs("Attribute");
  xmlAbs.addAttribute( "Active", "1" );
  xmlAbs.addAttribute( "AttributeType", "Scalar" );
  xmlAbs.addAttribute( "Center", "Node" );
  xmlAbs.addAttribute( "Name", "abs(psi)" );
  xmlAbs.addChild( xmlAbsData );


  // tell me where the actual ARG(PSI) data sits
  Teuchos::XMLObject xmlArgData("DataItem");
  xmlArgData.addAttribute( "NumberType", "Float" );
  xmlArgData.addAttribute( "Precision", "4" );
  os << "1 " << Nx+1 << " " << Nx+1;
  xmlArgData.addAttribute( "Dimensions", os.str() );
  os.str("");
  xmlArgData.addAttribute( "Format", "HDF" );
  xmlArgData.addContent( "myfile.h5:/MyGrid/arg(psi)" );

  Teuchos::XMLObject xmlArg("Attribute");
  xmlArg.addAttribute( "Center", "Node" );
  xmlArg.addAttribute( "AttributeType", "Scalar" );
  xmlArg.addAttribute( "Name", "arg(psi)" );
  xmlArg.addChild( xmlArgData );

  // put it all in GRID
  Teuchos::XMLObject xmlGrid("Grid");
  xmlGrid.addAttribute( "Name", "MyGrid" );
  xmlGrid.addChild( xmlTopology );
  xmlGrid.addChild( xmlGeometry );
  xmlGrid.addChild( xmlAbs );
  xmlGrid.addChild( xmlArg );

  // put grid in domain
  Teuchos::XMLObject xmlDomain("Domain");
  xmlDomain.addChild( xmlGrid );

  // out domain in Xdmf
  Teuchos::XMLObject xdmfContainer("Xdmf");
  xdmfContainer.addChild( xmlDomain );

  // open the file
  xdmfFile.open( filename.c_str() );
  // write the xml tree to a file
  xdmfFile << "<?xml version=\"1.0\" ?>" << std::endl
           << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [" << std::endl
           << "<!ENTITY HeavyData \"myfile.h5\">" << std::endl
           << "]>" << std::endl
           << std::endl;

  xdmfFile << xdmfContainer;
  // close the file
  xdmfFile.close();
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // write the HDF5 heavy data file
  // ---------------------------------------------------------------------------
  // create a vector with abs values
//   Epetra_MultiVector absPsi(StandardMap,1);
// 
// // std::cout << absPsi << std::endl;
// 
//   // fill absPsi
//   for (int i=0; i<Nx+1; i++) {
//       index[0] = i;
//       for (int j=0; j<Nx+1; j++) {
//           index[1] = j;
//           k = SGrid.i2k( index );
//           absPsi.ReplaceGlobalValue( k, 1, abs((*psi)[k]) );
//       }
//   }

//   EpetraExt::HDF5 myhdf5(comm);
//   myhdf5.Create("data/myfile.h5");
//   myhdf5.Write("MyGrid/abs(psi)", absPsi);
//   myhdf5.Write("MyGrid/arg(psi)", absPsi);
//   myhdf5.Close();
  // ---------------------------------------------------------------------------

}
// =============================================================================