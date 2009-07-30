#include "ioXdmf.h"

#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLParameterListWriter.hpp>
#include <EpetraExt_Utils.h> // for toString

#include <Epetra_MultiVector.h>
#include <EpetraExt_HDF5.h>

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// =============================================================================
// Constructor
IoXdmf::IoXdmf( StaggeredGrid &sGrid ):
  SGrid(sGrid) // make this &sGrid
{
}
// =============================================================================



// =============================================================================
// Destructor
IoXdmf::~IoXdmf()
{
}
// =============================================================================



// =============================================================================
void IoXdmf::read( const std::string           &fileName,
                   std::vector<double_complex> *psi,
                   Teuchos::ParameterList      *problemParams )
{
  std::cerr << "IoXdmf::read not yet implemented." << std::endl;
  exit(EXIT_FAILURE);
}
// =============================================================================



// =============================================================================
void IoXdmf::write( const std::string                 &fileName,
                    const std::vector<double_complex> &psi,
                    const Teuchos::ParameterList      &problemParams )
{

  int    Nx = SGrid.getNx(),
         k,
         index[2];
  double h  = SGrid.getH();
  std::string   str;
  std::ofstream xdmfFile;

  // create the HDF5 file name, take off the suffix and replace it by .h5
  int dotPos = fileName.rfind(".");
  std::string hdf5FileName = fileName;
  hdf5FileName.replace(dotPos,fileName.size()-dotPos ,".h5");
  // strip off the directory names to get the base name
  dotPos = fileName.rfind("/");
  std::string hdf5BaseName = hdf5FileName;
  hdf5BaseName.erase(0,dotPos+1);

  // open the file
  xdmfFile.open( fileName.c_str() );

  // write the XDMF header
  xdmfFile << "<?xml version=\"1.0\" ?>\n"
           << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
           << "<!ENTITY HeavyData \"" << hdf5BaseName << "\">\n"
           << "]>\n"
           << std::endl;

  // enclosing structure
  Teuchos::XMLObject xdmfContainer("Xdmf");

  // Append the problem parameters in XML form.
  // This needs to go *inside* the enclosing Xdmf tag, otherwise the file format
  // is not universally recognized.
  Teuchos::XMLObject xmlParameterList("");
  xmlParameterList = Teuchos::XMLParameterListWriter::XMLParameterListWriter()
                                                        .toXML( problemParams );
  xdmfContainer.addChild(xmlParameterList);

  // add the domain item
  Teuchos::XMLObject xmlDomain("Domain");
  xdmfContainer.addChild( xmlDomain );

  // add the grid item
  Teuchos::XMLObject xmlGrid("Grid");
  xmlGrid.addAttribute( "Name", "MyGrid" );
  xmlDomain.addChild( xmlGrid );

  // add topology
  Teuchos::XMLObject xmlTopology("Topology");
  xmlTopology.addAttribute( "TopologyType", "3DCORECTMESH" );
  str = "1 " + EpetraExt::toString(Nx+1) + " " + EpetraExt::toString(Nx+1);
  xmlTopology.addAttribute( "Dimensions", str );
  xmlGrid.addChild( xmlTopology );

  // merge the latter two into geometry
  Teuchos::XMLObject xmlGeometry("Geometry");
  xmlGeometry.addAttribute( "Type", "ORIGIN_DXDYDZ" );
  xmlGrid.addChild( xmlGeometry );

  // add actual geometry information to the geometry item
  // define origin
  Teuchos::XMLObject xmlOrigin("DataItem");
  xmlOrigin.addAttribute( "Name", "Origin" );
  xmlOrigin.addAttribute( "NumberType", "Float" );
  xmlOrigin.addAttribute( "Dimensions", "3" );
  xmlOrigin.addAttribute( "Format", "XML" );
  xmlOrigin.addContent( "0 0 0" );
  xmlGeometry.addChild( xmlOrigin );

  // define spacing
  Teuchos::XMLObject xmlSpacing("DataItem");
  xmlSpacing.addAttribute( "Name", "Spacing" );
  xmlSpacing.addAttribute( "NumberType", "Float" );
  xmlSpacing.addAttribute( "Dimensions", "3" );
  xmlSpacing.addAttribute( "Format", "XML" );
  str = "0 " + EpetraExt::toString(h) + " " + EpetraExt::toString(h);
  xmlSpacing.addContent( str );
  xmlGeometry.addChild( xmlSpacing );

  // tell me where the actual ABS(PSI) data sits
  Teuchos::XMLObject xmlAbs("Attribute");
  xmlAbs.addAttribute( "Active", "1" );
  xmlAbs.addAttribute( "AttributeType", "Scalar" );
  xmlAbs.addAttribute( "Center", "Node" );
  xmlAbs.addAttribute( "Name", "abs(psi)" );
  xmlGrid.addChild( xmlAbs );

  Teuchos::XMLObject xmlAbsData("DataItem");
  xmlAbsData.addAttribute( "NumberType", "Float" );
  xmlAbsData.addAttribute( "Precision", "4" );
  str = "1 " + EpetraExt::toString(Nx+1) + " " + EpetraExt::toString(Nx+1);
  xmlAbsData.addAttribute( "Dimensions", str );
  xmlAbsData.addAttribute( "Format", "HDF" );
  xmlAbsData.addContent( hdf5BaseName+":/abs/Values" );
  xmlAbs.addChild( xmlAbsData );


  // tell me where the actual ARG(PSI) data sits
  Teuchos::XMLObject xmlArg("Attribute");
  xmlArg.addAttribute( "Center", "Node" );
  xmlArg.addAttribute( "AttributeType", "Scalar" );
  xmlArg.addAttribute( "Name", "arg(psi)" );
  xmlGrid.addChild( xmlArg );

  Teuchos::XMLObject xmlArgData("DataItem");
  xmlArgData.addAttribute( "NumberType", "Float" );
  xmlArgData.addAttribute( "Precision", "4" );
  str = "1 " + EpetraExt::toString(Nx+1) + " " + EpetraExt::toString(Nx+1);
  xmlArgData.addAttribute( "Dimensions", str );
  xmlArgData.addAttribute( "Format", "HDF" );
  xmlArgData.addContent( hdf5BaseName+":/arg/Values" );
  xmlArg.addChild( xmlArgData );

  // write it all to the file
  xdmfFile << xdmfContainer;

  // close the file
  xdmfFile.close();
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // write the HDF5 heavy data file
  // ---------------------------------------------------------------------------
  // create a vector with abs values
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
//   int NumUnknowns = (SGrid.getNx()+1)*(SGrid.getNx()+1);
  Epetra_Map StandardMap( Nx+1,0,Comm);

  // fill absPsi and argPsi
  Epetra_MultiVector absPsi(StandardMap,Nx+1);
  Epetra_MultiVector argPsi(StandardMap,Nx+1);
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = SGrid.i2k( index );
          absPsi.ReplaceGlobalValue( i, j, abs(psi[k]) );
          argPsi.ReplaceGlobalValue( i, j, arg(psi[k]) );
      }
  }

  EpetraExt::HDF5 myhdf5( Comm );
  myhdf5.Create( hdf5FileName );
  myhdf5.Write( "abs", absPsi );
  myhdf5.Write( "arg", argPsi );
  myhdf5.Close();
  // ---------------------------------------------------------------------------

}
// =============================================================================