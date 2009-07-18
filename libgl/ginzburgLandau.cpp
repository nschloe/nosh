#include "ginzburgLandau.h"
#include <complex>
#include <iostream>


// =============================================================================
// Class constructor
GinzburgLandau::GinzburgLandau( int nx ):
Nx(nx),
d(2)
{
}
// =============================================================================


// =============================================================================
// Destructor
GinzburgLandau::~GinzburgLandau()
{
}
// =============================================================================


// // =============================================================================
// // calculate the 
// GinzburgLandau::computeF( psiReal, psiImag )
// {
// }
// // =============================================================================


// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
std::complex<double> GinzburgLandau::boundaryConditions( int eqnum,
                                                         double* psiReal,
                                                         double* psiImag,
                                                         PsiGrid::PsiGrid psiGrid,
                                                         AGrid::AGrid     aGrid )
{
  // Equations are ordered counter-clockwise, starting at the origin.
  static double x[2];

  // the preliminary result type  
  std::complex<double> res,
                       cUnit(0.0,1.0),
                       psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;
  double ARight, ALeft, AAbove, ABelow;               
  int i[2] = {0,0};
  int k;

  enum borderType
  {
     BOTTOMLEFT,
     BOTTOMRIGHT,
     TOPLEFT,
     TOPRIGHT,
     BOTTOM,
     TOP,
     LEFT,
     RIGHT
  };


  borderType bt;
  // define an order

  try{
      if (eqnum==0) {
          bt = BOTTOMLEFT;
          i[0] = 0; i[1] = 0;      
      } else if (eqnum==Nx) {
          bt = BOTTOMRIGHT;
          i[0] = Nx; i[1] = 0;
      } else if (eqnum==2*Nx) {
          bt = TOPRIGHT;
          i[0] = Nx; i[1] = Nx;
      } else if (eqnum==3*Nx) {
          bt = TOPLEFT;
          i[0] = 0; i[1] = Nx;
      } else if (eqnum<Nx) {
          bt = BOTTOM;
          i[0] = k; i[1] = 0;
      } else if (eqnum<2*Nx) {
          bt = RIGHT;
          i[0] = 0; i[1] = k-Nx;
      } else if (eqnum<3*Nx) {
          bt = TOP;
          i[0] = 3*Nx-k; i[1] = Nx;
      } else if (eqnum<4*Nx) {
          bt = LEFT;
          i[0] = 0; i[1] = 4*Nx-k;
      } else {
          throw 1;
      }
  }
  catch(int) {
      std::cout << "Illegal running index k=" << k << "in function boundaryConditions." << std::endl;
  }

  try {
      switch (bt) {
          case BOTTOMLEFT:
              k = psiGrid.i2k( i );
              psiK = ( psiReal[k], psiImag[k] );

              i[0] = i[0]+1;
              k = psiGrid.i2k( i );
              psiKRight = ( psiReal[k], psiImag[k] );
              ARight    = aGrid.getAxRight( i );

              i[1] = i[1]+1;
              k = psiGrid.i2k( i );
              psiKAbove = ( psiReal[k], psiImag[k] );
              AAbove    = aGrid.getAyAbove( i );

              res = ( - psiK * 2.0
                      + psiKRight * exp(-cUnit*ARight*h)
                      + psiKAbove * exp(-cUnit*AAbove*h) ) * cUnit / (sqrt(2)*h);
              break;

          case BOTTOMRIGHT:
              k = psiGrid.i2k( i );
              psiK = ( psiReal[k], psiImag[k] );

              i[0] = i[0]-1;
              k = psiGrid.i2k( i );
              psiKLeft = ( psiReal[k], psiImag[k] );
              ALeft    = aGrid.getAxLeft( i );

              i[1] = i[1]+1;
              k = psiGrid.i2k( i );
              psiKAbove = ( psiReal[k], psiImag[k] );
              AAbove    = aGrid.getAyAbove( i );

              res = ( - psiK * 2.0
                      + psiKLeft  * exp( cUnit*ALeft *h)
                      + psiKAbove * exp(-cUnit*AAbove*h) ) * cUnit / (sqrt(2)*h);
              break;

          case TOPRIGHT:
              k = psiGrid.i2k( i );
              psiK = ( psiReal[k], psiImag[k] );

              i[0] = i[0]-1;
              k = psiGrid.i2k( i );
              psiKLeft = ( psiReal[k], psiImag[k] );
              ALeft    = aGrid.getAxLeft( i );

              i[1] = i[1]-1;
              k = psiGrid.i2k( i );
              psiKBelow = ( psiReal[k], psiImag[k] );
              ABelow    = aGrid.getAyBelow( i );

              res = ( - psiK * 2.0
                      + psiKLeft  * exp( cUnit*ALeft *h)
                      + psiKBelow * exp( cUnit*ABelow*h) ) * cUnit / (sqrt(2)*h);
              break;

            
          case TOPLEFT:
              k = psiGrid.i2k( i );
              psiK = ( psiReal[k], psiImag[k] );

              i[0] = i[0]+1;
              k = psiGrid.i2k( i );
              psiKRight = ( psiReal[k], psiImag[k] );
              ARight    = aGrid.getAxRight( i );

              i[1] = i[1]-1;
              k = psiGrid.i2k( i );
              psiKBelow = ( psiReal[k], psiImag[k] );
              ABelow    = aGrid.getAyBelow( i );

              res = ( - psiK * 2.0
                      + psiKRight * exp( cUnit*ARight*h)
                      + psiKBelow * exp( cUnit*ABelow*h) ) * cUnit / (sqrt(2)*h);
              break;
            
          case BOTTOM:
              k = psiGrid.i2k( i );
              psiK = ( psiReal[k], psiImag[k] );

              i[1] = i[1]+1;
              k = psiGrid.i2k( i );
              psiKAbove = ( psiReal[k], psiImag[k] );
              AAbove    = aGrid.getAyAbove( i );

              res = ( - psiK
                      + psiKAbove * exp(-cUnit*AAbove*h) ) * cUnit/h;
              break;

          case RIGHT:
              k = psiGrid.i2k( i );
              psiK = ( psiReal[k], psiImag[k] );

              i[0] = i[0]-1;
              k = psiGrid.i2k( i );
              psiKLeft = ( psiReal[k], psiImag[k] );
              ALeft    = aGrid.getAxLeft( i );

              res = ( - psiK
                      + psiKLeft * exp(-cUnit*ALeft*h) ) * cUnit/h;
              break;

          case TOP:
              k = psiGrid.i2k( i );
              psiK = ( psiReal[k], psiImag[k] );

              i[1] = i[1]-1;
              k = psiGrid.i2k( i );
              psiKBelow = ( psiReal[k], psiImag[k] );
              ABelow    = aGrid.getAyBelow( i );

              res = ( - psiK
                      + psiKBelow * exp(-cUnit*ABelow*h) ) * cUnit/h;
              break;

          case LEFT:
              k = psiGrid.i2k( i );
              psiK = ( psiReal[k], psiImag[k] );

              i[0] = i[0]+1;
              k = psiGrid.i2k( i );
              psiKBelow = ( psiReal[k], psiImag[k] );
              ARight    = aGrid.getAxRight( i );

              res = ( - psiK
                      + psiKRight * exp(-cUnit*ARight*h) ) * cUnit/h;
              break;

          default:
              throw 1;
      }
  }
  catch (int) {
      std::cout << "Illegal enumerator in function boundaryConditions." << std::endl;
  }

  // return the result
  return res;
}
// =============================================================================



// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
std::complex<double> GinzburgLandau::interiorEquations( int eqnum,
                                                        std::complex<double>* psi,
                                                        PsiGrid::PsiGrid psiGrid,
                                                        AGrid::AGrid     aGrid )
{
  // Equations are ordered counter-clockwise, starting at the origin.
  static double x[2];

  std::complex<double> psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow,
                       res;
  const std::complex<double> cUnit(0.0,1.0);
  int i[2];
  int k;

  // get the corresponding running index for eqnum
  i[1] = eqnum%(Nx-1) + 1;
  i[2] = eqnum/(Nx-1) + 1;

  double ARight = aGrid.getAxRight( i );
  double ALeft  = aGrid.getAxRight( i );
  double ABelow = aGrid.getAyBelow( i );
  double AAbove = aGrid.getAyAbove( i );

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
          + exp(cUnit*ALeft *h) + exp(-cUnit*ARight*h)
          + exp(cUnit*ABelow*h) + exp(-cUnit*AAbove*h) ) / (h*h)
        + psiK * (1-abs(psiK)*abs(psiK));

  // return the result
  return res;
}
// =============================================================================