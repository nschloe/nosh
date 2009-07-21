#include "ginzburgLandau.h"
#include <iostream>
#include <vector>

// =============================================================================
// Class constructor
GinzburgLandau::GinzburgLandau( int nx,
                                double edgelength,
                                double h0 ):
h      ( edgelength/nx ),
psiGrid( PsiGrid::PsiGrid( nx,
                           edgelength ) ),
aGrid  ( AGrid::AGrid( nx,
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
// Defines a mapping of all the GL equations to a running index.
// The equations in the GL context are each associated with and
// uniquely identified by a particular node
// (around which the equation is centered, e.g., the center node with the five-
// point stencil).
void GinzburgLandau::getEquationType( int eqnum,
                                      equationType &eqType,
                                      int *i )
{
  int Nx = psiGrid.getNx();

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
std::complex<double> GinzburgLandau::computeGl( int eqnum,
                                                std::vector<std::complex<double> > psi )
{
  // Equations are ordered counter-clockwise, starting at the origin.

  // the preliminary result type  
  std::complex<double> res,
                       cUnit(0.0,1.0),
                       psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;
  double ARight, ALeft, AAbove, ABelow;
  int i[2];
  int k;
  equationType eqType;

  // get the equation type from the running index
  getEquationType( eqnum, eqType, i );

  switch (eqType) {
      case BOTTOMLEFT:
          k = psiGrid.i2k( i );
          psiK = psi[k];

          i[0] = i[0]+1;
          k = psiGrid.i2k( i );
          psiKRight = psi[k];
          ARight    = aGrid.getAxRight( i );
          i[0] = i[0]-1;

          i[1] = i[1]+1;
          k = psiGrid.i2k( i );
          psiKAbove = psi[k];
          AAbove    = aGrid.getAyAbove( i );
          i[1] = i[1]-1;

          res = ( - psiK      * 2.0
                  + psiKRight * exp(-cUnit*ARight*h)
                  + psiKAbove * exp(-cUnit*AAbove*h) ) * cUnit / (sqrt(2)*h);

          break;

      case BOTTOMRIGHT:
          k = psiGrid.i2k( i );
          psiK = psi[k];

          i[0] = i[0]-1;
          k = psiGrid.i2k( i );
          psiKLeft = psi[k];
          ALeft    = aGrid.getAxLeft( i );
          i[0] = i[0]+1;

          i[1] = i[1]+1;
          k = psiGrid.i2k( i );
          psiKAbove = psi[k];
          AAbove    = aGrid.getAyAbove( i );
          i[1] = i[1]-1;

          res = ( - psiK * 2.0
                  + psiKLeft  * exp( cUnit*ALeft *h)
                  + psiKAbove * exp(-cUnit*AAbove*h) ) * cUnit / (sqrt(2)*h);
          break;

      case TOPRIGHT:
          k = psiGrid.i2k( i );
          psiK = psi[k];

          i[0] = i[0]-1;
          k = psiGrid.i2k( i );
          psiKLeft = psi[k];
          ALeft    = aGrid.getAxLeft( i );
          i[0] = i[0]+1;

          i[1] = i[1]-1;
          k = psiGrid.i2k( i );
          psiKBelow = psi[k];
          ABelow    = aGrid.getAyBelow( i );
          i[1] = i[1]+1;

          res = ( - psiK * 2.0
                  + psiKLeft  * exp( cUnit*ALeft *h)
                  + psiKBelow * exp( cUnit*ABelow*h) ) * cUnit / (sqrt(2)*h);
          break;

      case TOPLEFT:
          k = psiGrid.i2k( i );
          psiK = psi[k];

          i[0] = i[0]+1;
          k = psiGrid.i2k( i );
          psiKRight = psi[k];
          ARight    = aGrid.getAxRight( i );
          i[0] = i[0]-1;

          i[1] = i[1]-1;
          k = psiGrid.i2k( i );
          psiKBelow = psi[k];
          ABelow    = aGrid.getAyBelow( i );
          i[1] = i[1]+1;

          res = ( - psiK * 2.0
                  + psiKRight * exp(-cUnit*ARight*h)
                  + psiKBelow * exp( cUnit*ABelow*h) ) * cUnit / (sqrt(2)*h);
          break;

      case BOTTOM:
          k = psiGrid.i2k( i );
          psiK = psi[k];

          i[1] = i[1]+1;
          k = psiGrid.i2k( i );
          psiKAbove = psi[k];
          AAbove    = aGrid.getAyAbove( i );
          i[1] = i[1]-1;

          res = ( - psiK
                  + psiKAbove * exp(-cUnit*AAbove*h) ) * cUnit/h;
          break;

      case RIGHT:
          k = psiGrid.i2k( i );
          psiK = psi[k];

          i[0] = i[0]-1;
          k = psiGrid.i2k( i );
          psiKLeft = psi[k];
          ALeft    = aGrid.getAxLeft( i );
          i[0] = i[0]+1;

          res = ( - psiK
                  + psiKLeft * exp(-cUnit*ALeft*h) ) * cUnit/h;
          break;

      case TOP:
          k = psiGrid.i2k( i );
          psiK = psi[k];

          i[1] = i[1]-1;
          k = psiGrid.i2k( i );
          psiKBelow = psi[k];
          ABelow    = aGrid.getAyBelow( i );
          i[1] = i[1]+1;

          res = ( - psiK
                  + psiKBelow * exp(-cUnit*ABelow*h) ) * cUnit/h;
          break;

      case LEFT:
          k = psiGrid.i2k( i );
          psiK = psi[k];

          i[0] = i[0]+1;
          k = psiGrid.i2k( i );
          psiKRight = psi[k];
          ARight    = aGrid.getAxRight( i );
          i[0] = i[0]-1;

          res = ( - psiK
                  + psiKRight * exp(-cUnit*ARight*h) ) * cUnit/h;
          break;

      case INTERIOR:
          ALeft  = aGrid.getAxLeft ( i );
          ARight = aGrid.getAxRight( i );
          ABelow = aGrid.getAyBelow( i );
          AAbove = aGrid.getAyAbove( i );

          k = psiGrid.i2k( i );
          psiK = psi[k];

          i[0] = i[0]+1;
          k = psiGrid.i2k( i );
          psiKRight = psi[k];
          i[0] = i[0]-1;

          i[0] = i[0]-1;
          k = psiGrid.i2k( i );
          psiKLeft = psi[k];
          i[0] = i[0]+1;

          i[1] = i[1]+1;
          k = psiGrid.i2k( i );
          psiKAbove = psi[k];
          i[1] = i[1]-1;

          i[1] = i[1]-1;
          k = psiGrid.i2k( i );
          psiKBelow = psi[k];
          i[1] = i[1]+1;

          res = ( -4.0*psiK
                  + psiKLeft*  exp(cUnit*ALeft *h) + psiKRight* exp(-cUnit*ARight*h)
                  + psiKBelow* exp(cUnit*ABelow*h) + psiKAbove* exp(-cUnit*AAbove*h) ) / (h*h)
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
void GinzburgLandau::getJacobianRow( int eqnum,
                                     std::vector<std::complex<double> > psi,
                                     std::vector<int>& columnIndicesPsi,
                                     std::vector<std::complex<double> >& valuesPsi,
                                     std::vector<int>& columnIndicesPsiConj,
                                     std::vector<std::complex<double> >& valuesPsiConj )
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
void GinzburgLandau::getJacobianRowSparsity( int eqnum,
                                             std::vector<int>& columnIndicesPsi,
                                             std::vector<int>& columnIndicesPsiConj )
{
  // create dummy arguments
  std::vector<std::complex<double> > psi;
  std::vector<std::complex<double> > valuesPsi,
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
void GinzburgLandau::computeJacobianRow( filltype ft,
                                         int eqnum,
                                         std::vector<std::complex<double> > psi,
                                         std::vector<int>& columnIndicesPsi,
                                         std::vector<std::complex<double> >& valuesPsi,
                                         std::vector<int>& columnIndicesPsiConj,
                                         std::vector<std::complex<double> >& valuesPsiConj )
{
  int i[2];
  int k, kLeft, kRight, kBelow, kAbove;
  int numEntriesPsi, numEntriesPsiConj;
  double ARight, ALeft, AAbove, ABelow;
  const std::complex<double> cUnit(0,1);
  std::complex<double> psiK, psiKLeft, psiKRight, psiKBelow, psiKAbove;
  equationType eqType;

  // get the equation type from the running index
  getEquationType( eqnum, eqType, i );

  switch (eqType) {
      case BOTTOMLEFT:
          k = psiGrid.i2k( i );
          if (ft==VALUES)
              psiK = psi[k];

          i[0] = i[0]+1;
          kRight    = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKRight = psi[kRight];
              ARight    = aGrid.getAxRight( i );
          }
          i[0] = i[0]-1;

          i[1] = i[1]+1;
          kAbove = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKAbove = psi[kAbove];
              AAbove    = aGrid.getAyAbove( i );
          }
          i[1] = i[1]-1;

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kRight;
          columnIndicesPsi[2] = kAbove;

          if (ft==VALUES) {
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 2.0                * cUnit / (sqrt(2)*h);
              valuesPsi[1] = exp(-cUnit*ARight*h) * cUnit / (sqrt(2)*h);
              valuesPsi[2] = exp(-cUnit*AAbove*h) * cUnit / (sqrt(2)*h);
          }

          break;

      case BOTTOMRIGHT:
          k = psiGrid.i2k( i );
          if (ft==VALUES)
              psiK = psi[k];

          i[0] = i[0]-1;
          kLeft = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKLeft = psi[kLeft];
              ALeft    = aGrid.getAxLeft( i );
          }
          i[0] = i[0]+1;

          i[1] = i[1]+1;
          kAbove = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKAbove = psi[kAbove];
              AAbove    = aGrid.getAyAbove( i );
          }
          i[1] = i[1]-1;

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kLeft;
          columnIndicesPsi[2] = kAbove;

          if (ft==VALUES) {
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 2.0                * cUnit / (sqrt(2)*h);
              valuesPsi[1] = exp( cUnit*ALeft *h) * cUnit / (sqrt(2)*h);
              valuesPsi[2] = exp(-cUnit*AAbove*h) * cUnit / (sqrt(2)*h);
          }

          break;

      case TOPRIGHT:
          k = psiGrid.i2k( i );
          if (ft==VALUES)
              psiK = psi[k];

          i[0] = i[0]-1;
          kLeft = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKLeft = psi[kLeft];
              ALeft    = aGrid.getAxLeft( i );
          }
          i[0] = i[0]+1;

          i[1] = i[1]-1;
          kBelow = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKBelow = psi[kBelow];
              ABelow    = aGrid.getAyBelow( i );
          }
          i[1] = i[1]+1;

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kLeft;
          columnIndicesPsi[2] = kBelow;

          if (ft==VALUES) {
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 2.0                * cUnit / (sqrt(2)*h);
              valuesPsi[1] = exp( cUnit*ALeft *h) * cUnit / (sqrt(2)*h);
              valuesPsi[2] = exp( cUnit*ABelow*h) * cUnit / (sqrt(2)*h);
          }

          break;

      case TOPLEFT:
          k = psiGrid.i2k( i );
          if (ft==VALUES)
              psiK = psi[k];

          i[0] = i[0]+1;
          kRight = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKRight = psi[kRight];
              ARight    = aGrid.getAxRight( i );
          }
          i[0] = i[0]-1;

          i[1] = i[1]-1;
          kBelow = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKBelow = psi[kBelow];
              ABelow    = aGrid.getAyBelow( i );
          }
          i[1] = i[1]+1;

          numEntriesPsi = 3;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kRight;
          columnIndicesPsi[2] = kBelow;

          if (ft==VALUES) {
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 2.0                * cUnit / (sqrt(2)*h);
              valuesPsi[1] = exp(-cUnit*ARight*h) * cUnit / (sqrt(2)*h);
              valuesPsi[2] = exp( cUnit*ABelow*h) * cUnit / (sqrt(2)*h);
          }

          break;

      case BOTTOM:
          k = psiGrid.i2k( i );
          if (ft==VALUES)
              psiK = psi[k];

          i[1] = i[1]+1;
          kAbove = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKAbove = psi[kAbove];
              AAbove    = aGrid.getAyAbove( i );
          }
          i[1] = i[1]-1;

          numEntriesPsi = 2;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kAbove;

          if (ft==VALUES) {
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 1.0                * cUnit/h;
              valuesPsi[1] = exp(-cUnit*AAbove*h) * cUnit/h;
          }

          break;

      case RIGHT:
          k = psiGrid.i2k( i );
          if (ft==VALUES)
              psiK = psi[k];

          i[0] = i[0]-1;
          kLeft = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKLeft = psi[kLeft];
              ALeft    = aGrid.getAxLeft( i );
          }
          i[0] = i[0]+1;

          numEntriesPsi = 2;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kLeft;

          if (ft==VALUES) {
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 1.0               * cUnit/h;
              valuesPsi[1] = exp(-cUnit*ALeft*h) * cUnit/h;
          }

          break;

      case TOP:
          k = psiGrid.i2k( i );
          if (ft==VALUES)
              psiK = psi[k];

          i[1] = i[1]-1;
          kBelow = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKBelow = psi[kBelow];
              ABelow    = aGrid.getAyBelow( i );
          }
          i[1] = i[1]+1;

          numEntriesPsi = 2;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kBelow;

          if (ft==VALUES) {
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 1.0                * cUnit/h;
              valuesPsi[1] = exp(-cUnit*ABelow*h) * cUnit/h;
          }

          break;

      case LEFT:
          k = psiGrid.i2k( i );
          if (ft==VALUES)
              psiK = psi[k];

          i[0] = i[0]+1;
          kRight = psiGrid.i2k( i );
          if (ft==VALUES) {
              psiKBelow = psi[k];
              ARight    = aGrid.getAxRight( i );
          }
          i[0] = i[0]-1;

          numEntriesPsi = 2;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kRight;

          if (ft==VALUES) {
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 1.0                * cUnit/h;
              valuesPsi[1] = exp(-cUnit*ARight*h) * cUnit/h;
          }
          break;

      case INTERIOR:
          ARight = aGrid.getAxRight( i );
          ALeft  = aGrid.getAxRight( i );
          ABelow = aGrid.getAyBelow( i );
          AAbove = aGrid.getAyAbove( i );

          k = psiGrid.i2k( i );
          if (ft==VALUES)
              psiK = psi[k];

          i[0] = i[0]+1;
          kRight = psiGrid.i2k( i );
          if (ft==VALUES)
              psiKRight = psi[k];
          i[0] = i[0]-1;

          i[0] = i[0]-1;
          kLeft = psiGrid.i2k( i );
          if (ft==VALUES)
              psiKLeft = psi[k];
          i[0] = i[0]+1;

          i[1] = i[1]+1;
          kAbove = psiGrid.i2k( i );
          if (ft==VALUES)
              psiKAbove = psi[k];
          i[1] = i[1]-1;

          i[1] = i[1]-1;
          kBelow = psiGrid.i2k( i );
          if (ft==VALUES)
              psiKBelow = psi[k];
          i[1] = i[1]+1;

          numEntriesPsi = 5;
          columnIndicesPsi.resize(numEntriesPsi);
          columnIndicesPsi[0] = k;
          columnIndicesPsi[1] = kLeft;
          columnIndicesPsi[2] = kRight;
          columnIndicesPsi[3] = kBelow;
          columnIndicesPsi[4] = kAbove;

          if (ft==VALUES) {
              valuesPsi.resize(numEntriesPsi);
              valuesPsi[0] = - 4.0                / (h*h)
                          + (1 - 2.0*abs(psiK)*abs(psiK));
              valuesPsi[1] = exp( cUnit*ALeft *h) / (h*h);
              valuesPsi[2] = exp(-cUnit*ARight*h) / (h*h);
              valuesPsi[3] = exp( cUnit*ABelow*h) / (h*h);
              valuesPsi[4] = exp(-cUnit*AAbove*h) / (h*h);
          }

          numEntriesPsiConj = 1;
          columnIndicesPsiConj.resize(numEntriesPsiConj);
          columnIndicesPsiConj[0] = k;

          if (ft==VALUES) {
              valuesPsiConj.resize(numEntriesPsiConj);
              valuesPsiConj[0] = -psiK*psiK;
          }
          break;

      default:
          std::cerr << "getJacobianRow:" << std::endl
                    << "    Illegal equationType \"" << eqType << "\". Abort." << std::endl;
          exit(EXIT_FAILURE);
  }
}
// =============================================================================