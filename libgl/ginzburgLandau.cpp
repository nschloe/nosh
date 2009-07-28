#include "ginzburgLandau.h"
#include <iostream>

#include <fstream> // for the VTK writer
#include <string>

#include <vector>

#include <Teuchos_XMLObject.hpp> // for the XMLized VTK format

#include <EpetraExt_HDF5.h> // for output in XDMF format
// #include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>


// #include <Teuchos_ParameterList.hpp>

#include <Teuchos_XMLParameterListWriter.hpp>
#include <Teuchos_XMLParameterListReader.hpp>


#include <EpetraExt_Utils.h> // for the toString function

#include <Epetra_SerialComm.h>
#include <EpetraExt_HDF5.h>

#include <Teuchos_FileInputSource.hpp> // to read XML data from files

// complex unit
const double_complex I(0,1);

// =============================================================================
// Class constructor
GinzburgLandau::GinzburgLandau( int    nx,
                                double edgelength,
                                double h0          ):
  sGrid( StaggeredGrid::StaggeredGrid( nx,
                                       edgelength,
                                       h0          ) )
{
}
// =============================================================================


// =============================================================================
// Destructor
GinzburgLandau::~GinzburgLandau()
{
}
// =============================================================================


// =============================================================================
StaggeredGrid::StaggeredGrid* GinzburgLandau::getStaggeredGrid()
{
    return &sGrid;
}
// =============================================================================


// =============================================================================
// Defines a mapping of all the GL equations to a running index.
// The equations in the GL context are each associated with and
// uniquely identified by a particular node
// (around which the equation is centered, e.g., the center node with the five-
// point stencil).
void GinzburgLandau::getEquationType( const int    eqnum,
                                      equationType &eqType,
                                      int          *i       )
{
  int Nx = sGrid.getNx();

  if (eqnum==0) {
      eqType = BOTTOMLEFT;
      i[0] = 0;
      i[1] = 0;
  } else if (eqnum==Nx) {
      eqType = BOTTOMRIGHT;
      i[0] = Nx;
      i[1] = 0;
  } else if (eqnum==2*Nx) {
      eqType = TOPRIGHT;
      i[0] = Nx;
      i[1] = Nx;
  } else if (eqnum==3*Nx) {
      eqType = TOPLEFT;
      i[0] = 0;
      i[1] = Nx;
  } else if (eqnum<Nx) {
      eqType = BOTTOM;
      i[0] = eqnum;
      i[1] = 0;
  } else if (eqnum<2*Nx) {
      eqType = RIGHT;
      i[0] = Nx;
      i[1] = eqnum-Nx;
  } else if (eqnum<3*Nx) {
      eqType = TOP;
      i[0] = 3*Nx-eqnum;
      i[1] = Nx;
  } else if (eqnum<4*Nx) {
      eqType = LEFT;
      i[0] = 0;
      i[1] = 4*Nx-eqnum;
  } else if (eqnum<(Nx+1)*(Nx+1)) {
      eqType = INTERIOR;
      i[0] = (eqnum-4*Nx)%(Nx-1) + 1;
      i[1] = (eqnum-4*Nx)/(Nx-1) + 1;
  } else {
      std::cerr << "getEquationType:" << std::endl
                << "    Illegal running index eqnum=" << eqnum << "." << std::endl;
      exit(EXIT_FAILURE);
  }
}
// =============================================================================


// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
double_complex GinzburgLandau::computeGl( const int                         eqnum,
                                          const std::vector<double_complex> &psi    )
{
  // Equations are ordered counter-clockwise, starting at the origin.

  // the preliminary result type
  double_complex res,
                 psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;
  double ARight, ALeft, AAbove, ABelow,
         h = sGrid.getH();
  int i[2];
  int k;
  equationType eqType;

  // get the equation type from the running index
  getEquationType( eqnum, eqType, i );

  switch (eqType) {
      case BOTTOMLEFT:
          k = sGrid.i2k( i );
          psiK = psi[k];

          i[0] += 1;
          k = sGrid.i2k( i );
          psiKRight = psi[k];
          i[0] -= 1;

          i[1] += 1;
          k = sGrid.i2k( i );
          psiKAbove = psi[k];
          i[1] -= 1;

          AAbove = sGrid.getAyAbove( i );
          ARight = sGrid.getAxRight( i );

          res = ( - psiK      * 2.0
                  + psiKRight * exp(-I*ARight*h)
                  + psiKAbove * exp(-I*AAbove*h) ) * I/(sqrt(2)*h);

          break;

      case BOTTOMRIGHT:
          k = sGrid.i2k( i );
          psiK = psi[k];

          i[0] -= 1;
          k = sGrid.i2k( i );
          psiKLeft = psi[k];
          i[0] += 1;

          i[1] += 1;
          k = sGrid.i2k( i );
          psiKAbove = psi[k];
          i[1] -= 1;

          ALeft  = sGrid.getAxLeft ( i );
          AAbove = sGrid.getAyAbove( i );

          res = ( - psiK      * 2.0
                  + psiKLeft  * exp( I*ALeft *h)
                  + psiKAbove * exp(-I*AAbove*h) ) * I/(sqrt(2)*h);
          break;

      case TOPRIGHT:
          k = sGrid.i2k( i );
          psiK = psi[k];

          i[0] -= 1;
          k = sGrid.i2k( i );
          psiKLeft = psi[k];
          i[0] += 1;

          i[1] -= 1;
          k = sGrid.i2k( i );
          psiKBelow = psi[k];
          i[1] += 1;

          ALeft  = sGrid.getAxLeft( i );
          ABelow = sGrid.getAyBelow( i );

          res = ( - psiK      * 2.0
                  + psiKLeft  * exp( I*ALeft *h)
                  + psiKBelow * exp( I*ABelow*h) ) * I/(sqrt(2)*h);
          break;

      case TOPLEFT:
          k = sGrid.i2k( i );
          psiK = psi[k];

          i[0] += 1;
          k = sGrid.i2k( i );
          psiKRight = psi[k];
          i[0] -= 1;

          i[1] -= 1;
          k = sGrid.i2k( i );
          psiKBelow = psi[k];
          i[1] += 1;

          ARight = sGrid.getAxRight( i );
          ABelow = sGrid.getAyBelow( i );

          res = ( - psiK      * 2.0
                  + psiKRight * exp(-I*ARight*h)
                  + psiKBelow * exp( I*ABelow*h) ) * I/(sqrt(2)*h);
          break;

      case BOTTOM:
          k = sGrid.i2k( i );
          psiK = psi[k];

          i[1] += 1;
          k = sGrid.i2k( i );
          psiKAbove = psi[k];
          i[1] -= 1;

          AAbove = sGrid.getAyAbove( i );

          res = ( - psiK
                  + psiKAbove * exp(-I*AAbove*h) ) * I/h;
          break;

      case RIGHT:
          k = sGrid.i2k( i );
          psiK = psi[k];

          i[0] -= 1;
          k = sGrid.i2k( i );
          psiKLeft = psi[k];
          i[0] += 1;

          ALeft = sGrid.getAxLeft( i );

          res = ( - psiK
                  + psiKLeft * exp( I*ALeft*h) ) * I/h;
          break;

      case TOP:
          k = sGrid.i2k( i );
          psiK = psi[k];

          i[1] -= 1;
          k = sGrid.i2k( i );
          psiKBelow = psi[k];
          i[1] += 1;

          ABelow = sGrid.getAyBelow( i );

          res = ( - psiK
                  + psiKBelow * exp( I*ABelow*h) ) * I/h;
          break;

      case LEFT:
          k = sGrid.i2k( i );
          psiK = psi[k];

          i[0] += 1;
          k = sGrid.i2k( i );
          psiKRight = psi[k];
          i[0] -= 1;

          ARight = sGrid.getAxRight( i );

          res = ( - psiK
                  + psiKRight * exp(-I*ARight*h) ) * I/h;
          break;

      case INTERIOR:
          k = sGrid.i2k( i );
          psiK = psi[k];

          i[0] += 1;
          k = sGrid.i2k( i );
          psiKRight = psi[k];
          i[0] -= 1;

          i[0] -= 1;
          k = sGrid.i2k( i );
          psiKLeft = psi[k];
          i[0] += 1;

          i[1] += 1;
          k = sGrid.i2k( i );
          psiKAbove = psi[k];
          i[1] -= 1;

          i[1] -= 1;
          k = sGrid.i2k( i );
          psiKBelow = psi[k];
          i[1] += 1;

          ALeft  = sGrid.getAxLeft ( i );
          ARight = sGrid.getAxRight( i );
          ABelow = sGrid.getAyBelow( i );
          AAbove = sGrid.getAyAbove( i );

          res = (   psiK*      (-4.0)
                  + psiKLeft*  exp( I*ALeft *h) + psiKRight* exp(-I*ARight*h)
                  + psiKBelow* exp( I*ABelow*h) + psiKAbove* exp(-I*AAbove*h) ) / (h*h)
                + psiK * (1-abs(psiK)*abs(psiK));
          break;

      default:
          std::cerr << "computeGl:" << std::endl
                    << "    Illegal equationType \"" << eqType << "\". Abort." << std::endl;
          exit(EXIT_FAILURE);
  }

  // return the result
  return res;
}
// =============================================================================



// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
void GinzburgLandau::getJacobianRow( const int                         eqnum,
                                     const std::vector<double_complex> &psi,
                                     std::vector<int>                  &columnIndicesPsi,
                                     std::vector<double_complex>       &valuesPsi,
                                     std::vector<int>                  &columnIndicesPsiConj,
                                     std::vector<double_complex>       &valuesPsiConj )
{
  computeJacobianRow( VALUES,
                      eqnum,
                      psi,
                      columnIndicesPsi,
                      valuesPsi,
                      columnIndicesPsiConj,
                      valuesPsiConj );

  return;
}
// =============================================================================



// =============================================================================
void GinzburgLandau::getJacobianRowSparsity( const int        eqnum,
                                             std::vector<int> &columnIndicesPsi,
                                             std::vector<int> &columnIndicesPsiConj )
{
  // create dummy arguments
  std::vector<double_complex> psi;
  std::vector<double_complex> valuesPsi,
                              valuesPsiConj;

  computeJacobianRow( SPARSITY,
                      eqnum,
                      psi,
                      columnIndicesPsi,
                      valuesPsi,
                      columnIndicesPsiConj,
                      valuesPsiConj );

  return;
}
// =============================================================================



// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
void GinzburgLandau::computeJacobianRow( const filltype                    ft,
                                         const int                         eqnum,
                                         const std::vector<double_complex> &psi,
                                         std::vector<int>                  &columnIndicesPsi,
                                         std::vector<double_complex>       &valuesPsi,
                                         std::vector<int>                  &columnIndicesPsiConj,
                                         std::vector<double_complex>       &valuesPsiConj )
{
  int          i[2],
               k, kLeft, kRight, kBelow, kAbove,
               numEntriesPsi, numEntriesPsiConj;
  double       ARight, ALeft, AAbove, ABelow,
               h = sGrid.getH();
  equationType eqType;

  // get the equation type from the running index
  getEquationType( eqnum, eqType, i );

  switch (eqType) {
      case BOTTOMLEFT:
          k = sGrid.i2k( i );

          i[0] += 1;
          kRight = sGrid.i2k( i );
          i[0] -= 1;

          i[1] += 1;
          kAbove = sGrid.i2k( i );
          i[1] -= 1;

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kRight;
          columnIndicesPsi[2] = kAbove;

          if (ft==VALUES) {
              ARight = sGrid.getAxRight( i );
              AAbove = sGrid.getAyAbove( i );
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = -2.0             * I/(sqrt(2)*h);
              valuesPsi[1] = exp(-I*ARight*h) * I/(sqrt(2)*h);
              valuesPsi[2] = exp(-I*AAbove*h) * I/(sqrt(2)*h);
          }

          break;

      case BOTTOMRIGHT:
          k = sGrid.i2k( i );

          i[0] -= 1;
          kLeft = sGrid.i2k( i );
          i[0] += 1;

          i[1] += 1;
          kAbove = sGrid.i2k( i );
          i[1] -= 1;

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kLeft;
          columnIndicesPsi[2] = kAbove;

          if (ft==VALUES) {
              ALeft    = sGrid.getAxLeft ( i );
              AAbove   = sGrid.getAyAbove( i );
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = -2.0             * I/(sqrt(2)*h);
              valuesPsi[1] = exp( I*ALeft *h) * I/(sqrt(2)*h);
              valuesPsi[2] = exp(-I*AAbove*h) * I/(sqrt(2)*h);
          }

          break;

      case TOPRIGHT:
          k = sGrid.i2k( i );

          i[0] -= 1;
          kLeft = sGrid.i2k( i );
          i[0] += 1;

          i[1] -= 1; 
          kBelow = sGrid.i2k( i );
          i[1] += 1;

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kLeft;
          columnIndicesPsi[2] = kBelow;

          if (ft==VALUES) {
              ALeft    = sGrid.getAxLeft ( i );
              ABelow   = sGrid.getAyBelow( i );
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = -2.0             * I/(sqrt(2)*h);
              valuesPsi[1] = exp( I*ALeft *h) * I/(sqrt(2)*h);
              valuesPsi[2] = exp( I*ABelow*h) * I/(sqrt(2)*h);
          }

          break;

      case TOPLEFT:
          k = sGrid.i2k( i );

          i[0] += 1;
          kRight = sGrid.i2k( i );
          i[0] -= 1;

          i[1] -= 1;
          kBelow = sGrid.i2k( i );
          i[1] += 1;

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kRight;
          columnIndicesPsi[2] = kBelow;

          if (ft==VALUES) {
              ARight    = sGrid.getAxRight( i );
              ABelow    = sGrid.getAyBelow( i );
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = -2.0             * I/(sqrt(2)*h);
              valuesPsi[1] = exp(-I*ARight*h) * I/(sqrt(2)*h);
              valuesPsi[2] = exp( I*ABelow*h) * I/(sqrt(2)*h);
          }

          break;

      case BOTTOM:
          k = sGrid.i2k( i );

          i[1] += 1;
          kAbove = sGrid.i2k( i );
          i[1] -= 1;

          numEntriesPsi = 2;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kAbove;

          if (ft==VALUES) {
              AAbove = sGrid.getAyAbove( i );
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = -1.0             * I/h;
              valuesPsi[1] = exp(-I*AAbove*h) * I/h;
          }

          break;

      case RIGHT:
          k = sGrid.i2k( i );

          i[0] -= 1;
          kLeft = sGrid.i2k( i );
          i[0] += 1;

          numEntriesPsi = 2;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kLeft;

          if (ft==VALUES) {
              ALeft = sGrid.getAxLeft( i );
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = -1.0            * I/h;
              valuesPsi[1] = exp( I*ALeft*h) * I/h;
          }

          break;

      case TOP:
          k = sGrid.i2k( i );

          i[1] -= 1;
          kBelow = sGrid.i2k( i );
          i[1] += 1;

          numEntriesPsi = 2;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kBelow;

          if (ft==VALUES) {
              ABelow = sGrid.getAyBelow( i );
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = -1.0             * I/h;
              valuesPsi[1] = exp( I*ABelow*h) * I/h;
          }

          break;

      case LEFT:
          k = sGrid.i2k( i );

          i[0] += 1;
          kRight = sGrid.i2k( i );
          i[0] -= 1;

          numEntriesPsi = 2;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kRight;

          if (ft==VALUES) {
              ARight = sGrid.getAxRight( i );
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = -1.0             * I/h;
              valuesPsi[1] = exp(-I*ARight*h) * I/h;
          }
          break;

      case INTERIOR:
          k = sGrid.i2k( i );

          i[0] += 1;
          kRight = sGrid.i2k( i );
          i[0] -= 1;

          i[0] -= 1;
          kLeft = sGrid.i2k( i );
          i[0] += 1;

          i[1] += 1;
          kAbove = sGrid.i2k( i );
          i[1] -= 1;

          i[1] -= 1;
          kBelow = sGrid.i2k( i );
          i[1] += 1;

          numEntriesPsi = 5;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kLeft;
          columnIndicesPsi[2] = kRight;
          columnIndicesPsi[3] = kBelow;
          columnIndicesPsi[4] = kAbove;

          if (ft==VALUES) {
              ARight = sGrid.getAxRight( i );
              ALeft  = sGrid.getAxLeft ( i );
              ABelow = sGrid.getAyBelow( i );
              AAbove = sGrid.getAyAbove( i );

              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 4.0            / (h*h)
                             + (1 - 2.0*abs(psi[k])*abs(psi[k]));
              valuesPsi[1] = exp( I*ALeft *h) / (h*h);
              valuesPsi[2] = exp(-I*ARight*h) / (h*h);
              valuesPsi[3] = exp( I*ABelow*h) / (h*h);
              valuesPsi[4] = exp(-I*AAbove*h) / (h*h);
          }

          numEntriesPsiConj = 1;
          columnIndicesPsiConj.resize(numEntriesPsiConj);
          columnIndicesPsiConj[0] = k;

          if (ft==VALUES) {
              valuesPsiConj.resize(numEntriesPsiConj);
              valuesPsiConj[0] = -psi[k]*psi[k];
          }
          break;

      default:
          std::cerr << "getJacobianRow:" << std::endl
                    << "    Illegal equationType \"" << eqType << "\". Abort." << std::endl;
          exit(EXIT_FAILURE);
  }
}
// =============================================================================


// =============================================================================
// calculate the free energy of a state
double GinzburgLandau::freeEnergy( const std::vector<double_complex> &psi )
{
  double                  energy = 0.0,
                          h      = sGrid.getH();
  StaggeredGrid::nodeType nt;

  for (unsigned int k=0; k<psi.size(); k++) {
      nt = sGrid.k2nodeType(k);
      if (nt==StaggeredGrid::CORNER)
          energy +=       h*h * pow(abs(psi[k]),4);
      else if (nt==StaggeredGrid::EDGE)
          energy += 0.5*  h*h * pow(abs(psi[k]),4);
      else if (nt==StaggeredGrid::INTERIOR)
          energy += 0.25* h*h * pow(abs(psi[k]),4);
      else
          std::cerr << "GinzburgLandau::freeEnergy" << std::endl
                    << "    Illegal node type " << nt << ". Abort." << std::endl;
  }

  return energy;
}
// =============================================================================


// =============================================================================
// calculate the free energy of a state
void GinzburgLandau::psiToLegacyVtkFile( const std::vector<double_complex> &psi,
                                         const std::string                 &filename )
{
  int           Nx = sGrid.getNx(),
                k,
                index[2];
  double        h  = sGrid.getH();
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
          << "POINT_DATA " << sGrid.getNumComplexUnknowns() << "\n";

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
          k = sGrid.i2k( index );
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
          k = sGrid.i2k( index );
          vtkfile << arg(psi[k]) << "\n";
      }
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // close the file
  vtkfile.close();
}
// =============================================================================



// =============================================================================
void GinzburgLandau::psiToVtkFile( const std::vector<double_complex> &psi,
                                   const Teuchos::ParameterList      &problemParams,
                                   const std::string                 &filename )
{
  int    Nx = sGrid.getNx(),
         k,
         index[2];
  double h  = sGrid.getH();

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
          k = sGrid.i2k( index );
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
          k = sGrid.i2k( index );
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
// parses and XML-style VTK file and returns psi as well as the problem
// parameters stored in the file
void GinzburgLandau::vtkFileToPsi( const std::string           &filename )
//                                    std::vector<double_complex> *psi,
//                                    Teuchos::ParameterList      *problemParams
{

// pass a possible 
//<?xml version="1.0"?>
// at the beginning of the file

std::cout << "11" << std::endl;
  Teuchos::FileInputSource xmlFile(filename);

std::cout << filename << std::endl;

  // extract the object from the filename
  Teuchos::XMLObject xmlFileObject = xmlFile.getObject();

  // plot the contents
  std::cout << xmlFileObject << std::endl;

  Teuchos::ParameterList plist;

std::cout << "22" << std::endl;
  // loop over the children
  const Teuchos::XMLObject* currentNode = &xmlFileObject;

  // find and read the parameter list
  currentNode = xmlBeagle ( &xmlFileObject, "ParameterList" );
  plist = Teuchos::XMLParameterListReader().toParameterList( *currentNode );

  std::cout << " plist: " << plist << std::endl;

std::cout << "33" << std::endl;

}
// =============================================================================



// =============================================================================
// Inside an XML object, this function looks for a specific tag and returns
// a pointer to it.
const Teuchos::XMLObject* GinzburgLandau::xmlBeagle ( const Teuchos::XMLObject *xmlObj,
                                                      const std::string        tag     )
{
  const Teuchos::XMLObject* xmlOut=NULL;

  if ( !xmlObj->getTag().compare(tag) ) // strings are equal
      return xmlObj;
  else
      for (int k=0; k<xmlObj->numChildren(); k++) {
          xmlOut = GinzburgLandau::xmlBeagle ( &(xmlObj->getChild(k)), tag );
          if (xmlOut) break; // not the null pointer => return
      }

  return xmlOut;
}
// =============================================================================



// =============================================================================
void GinzburgLandau::psiToXdmfFile( const std::vector<double_complex> &psi,
                                    const std::string                 &filename,
                                    const Epetra_Map                  &StandardMap,
                                    const Epetra_Comm                 &comm )
{

// TODO: use toString from EpetraExt instead of those nasty ofstreams!

  int    Nx = sGrid.getNx(),
         k,
         index[2];
  double h  = sGrid.getH();
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
  Epetra_MultiVector absPsi(StandardMap,1);

// std::cout << absPsi << std::endl;

  // fill absPsi
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = sGrid.i2k( index );
          absPsi.ReplaceGlobalValue( k, 1, abs(psi[k]) );
      }
  }

//   EpetraExt::HDF5 myhdf5(comm);
//   myhdf5.Create("data/myfile.h5");
//   myhdf5.Write("MyGrid/abs(psi)", absPsi);
//   myhdf5.Write("MyGrid/arg(psi)", absPsi);
//   myhdf5.Close();
  // ---------------------------------------------------------------------------

}
// =============================================================================