#include "ginzburgLandau.h"
#include <iostream>

#include <fstream> // for the VTK writer
#include <string>

#include <vector>

#include <Teuchos_XMLObject.hpp> // for the XMLized VTK format

// complex unit
const double_complex I(0,1);

// =============================================================================
// Class constructor
GinzburgLandau::GinzburgLandau( int    nx,
                                double edgelength,
                                double h0 ):
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
// calculate the free energy of a state
void GinzburgLandau::psiToVtkFile( const std::vector<double_complex> &psi,
                                   const std::string                 &filename )
{
  int    Nx = sGrid.getNx(),
         k,
         index[2];
  double h  = sGrid.getH();
  char*  str;

  std::ofstream      vtkfile;
  Teuchos::XMLObject vtuxml("VTKFile");

  // first build the XML structure
  Teuchos::XMLObject xmlDataArray("DataArray");
  xmlDataArray.addAttribute( "type", "Float32" );
  xmlDataArray.addAttribute( "name", "abs(psi)" );
  xmlDataArray.addAttribute( "format", "ascii" );
  for (int i=0; i<Nx+1; i++) {
      index[0] = i;
      for (int j=0; j<Nx+1; j++) {
          index[1] = j;
          k = sGrid.i2k( index );
          sprintf( str,"%f ",abs(psi[k]));
          xmlDataArray.addContent( str );
      }
  }

  Teuchos::XMLObject xmlPointData("PointData");
  xmlPointData.addAttribute( "Scalars", "abs(psi)" );
  xmlPointData.addChild(xmlDataArray);

  Teuchos::XMLObject xmlPiece("Piece");
  sprintf( str, "0 %d 0 %d 0 0", Nx, Nx );
  xmlPiece.addAttribute( "Extent", str );
  xmlPiece.addChild(xmlPointData);

  Teuchos::XMLObject xmlImageData("ImageData");
  sprintf( str, "0 %d 0 %d 0 0", Nx, Nx );
  xmlImageData.addAttribute( "WholeExtent", str );
  xmlImageData.addAttribute( "Origin", "0 0 0" );
  sprintf( str, "%f %f 0", h, h );
  xmlImageData.addAttribute( "Spacing", str );
  xmlImageData.addChild(xmlPiece);

  vtuxml.addAttribute( "type", "ImageData" );
  vtuxml.addAttribute( "version", "0.1" );
  vtuxml.addAttribute( "byte_order", "LittleEndian" );
  vtuxml.addChild(xmlImageData);

  // open the file
  vtkfile.open( filename.c_str() );
  // write the xml tree to a file
  vtkfile << "<?xml version=\"1.0\"?>" << std::endl;
  vtkfile << vtuxml;
  // close the file
  vtkfile.close();

std::cout << "done8" << std::endl;
}
// =============================================================================
