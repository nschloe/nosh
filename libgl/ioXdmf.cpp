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

// TODO: use toString from EpetraExt instead of those nasty ofstreams!

  int    Nx = SGrid.getNx(),
         k,
         index[2];
  double h  = SGrid.getH();
  std::ostringstream os;
  std::ofstream      xdmfFile;

  // create the HDF5 file name, take off the suffix and replace it by .h5
  int dotPos = fileName.rfind(".");
  std::string hdf5FileName = fileName;
  hdf5FileName.replace(dotPos,fileName.size()-dotPos ,".h5");
  // strip off the directory names to get the base name
  dotPos = fileName.rfind("/");
  std::string hdf5BaseName = hdf5FileName;
  hdf5BaseName.erase(0,dotPos+1);

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
  xmlAbsData.addContent( hdf5BaseName+":/abs/Values" );

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
  xmlArgData.addContent( hdf5BaseName+":/arg/Values" );

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
  xdmfFile.open( fileName.c_str() );
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