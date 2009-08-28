#include "ginzburgLandau.h"
#include "glException.h"

#include <iostream>

#include <fstream> // for the VTK writer
#include <string>

#include <vector>

#include <Teuchos_XMLObject.hpp> // for the XMLized VTK format

#include <EpetraExt_HDF5.h> // for output in XDMF format
// #include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>


#include <Teuchos_XMLParameterListWriter.hpp>
#include <Teuchos_XMLParameterListReader.hpp>


#include <EpetraExt_Utils.h> // for the toString function

#include <Epetra_SerialComm.h>
#include <EpetraExt_HDF5.h>

#include <Teuchos_FileInputSource.hpp> // to read XML data from files

// complex unit
const double_complex I ( 0,1 );

// =============================================================================
// Class constructor
GinzburgLandau::GinzburgLandau ( int    nx,
                                 double edgelength,
                                 double h0 ) :
    sGrid ( StaggeredGrid::StaggeredGrid ( nx,
                                           edgelength,
                                           h0 ) )
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
void GinzburgLandau::getEquationType ( const int    eqnum,
                                       equationType &eqType,
                                       int          *i )
{
  int Nx = sGrid.getNx();

  if ( eqnum==0 )
    {
      eqType = BOTTOMLEFT;
      i[0] = 0;
      i[1] = 0;
    }
  else if ( eqnum==Nx )
    {
      eqType = BOTTOMRIGHT;
      i[0] = Nx;
      i[1] = 0;
    }
  else if ( eqnum==2*Nx )
    {
      eqType = TOPRIGHT;
      i[0] = Nx;
      i[1] = Nx;
    }
  else if ( eqnum==3*Nx )
    {
      eqType = TOPLEFT;
      i[0] = 0;
      i[1] = Nx;
    }
  else if ( eqnum<Nx )
    {
      eqType = BOTTOM;
      i[0] = eqnum;
      i[1] = 0;
    }
  else if ( eqnum<2*Nx )
    {
      eqType = RIGHT;
      i[0] = Nx;
      i[1] = eqnum-Nx;
    }
  else if ( eqnum<3*Nx )
    {
      eqType = TOP;
      i[0] = 3*Nx-eqnum;
      i[1] = Nx;
    }
  else if ( eqnum<4*Nx )
    {
      eqType = LEFT;
      i[0] = 0;
      i[1] = 4*Nx-eqnum;
    }
  else if ( eqnum< ( Nx+1 ) * ( Nx+1 ) )
    {
      eqType = INTERIOR;
      i[0] = ( eqnum-4*Nx ) % ( Nx-1 ) + 1;
      i[1] = ( eqnum-4*Nx ) / ( Nx-1 ) + 1;
    }
  else
    {
      throw glException ( "GinzburgLandau::getEquationType",
                          "Illegal running index eqnum="
                          + EpetraExt::toString ( eqnum ) + "." );
    }
}
// =============================================================================


// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
double_complex GinzburgLandau::computeGl ( const int                         eqnum,
    const std::vector<double_complex> &psi )
{
  // Equations are ordered counter-clockwise, starting at the origin.

  // the preliminary result type
  double_complex res,
  psiK, psiKRight, psiKLeft, psiKAbove, psiKBelow;
  double ARight, ALeft, AAbove, ABelow;
  double h = sGrid.getH();
  int i[2];
  equationType eqType;

  // get the equation type from the running index
  getEquationType ( eqnum, eqType, i );

  // set psiK here, it is needed everywhere anyway
  psiK = psi[ sGrid.i2k ( i ) ];

  switch ( eqType )
    {
    case BOTTOMLEFT:
//           // -------------------------------------------------------------------
//           psiKRight = psi[ sGrid.getKRight(i) ];
//           psiKAbove = psi[ sGrid.getKAbove(i) ];
//
//           ARight = sGrid.getAxRight( i );
//           AAbove = sGrid.getAyAbove( i );
//
//           res = ( - psiK      * 2.0
//                   + psiKRight * exp(-I*ARight*h)
//                   + psiKAbove * exp(-I*AAbove*h) ) * I/(sqrt(2)*h);
//           // -------------------------------------------------------------------

      // -------------------------------------------------------------------
      psiKRight = psi[ sGrid.getKRight ( i ) ];
      psiKAbove = psi[ sGrid.getKAbove ( i ) ];

      ARight = sGrid.getAxRight ( i );
      AAbove = sGrid.getAyAbove ( i );

//       res = ( psiK* ( -2.0 )
//               + psiKRight* exp ( -I*ARight*h )
//               + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
//           + psiK * (1-norm(psiK));
      res = (   psiK     * (-2.0) ) / (h*h)
          + psiKRight * (-30.0)
          + psiK * (1-norm(psiK));
      // -------------------------------------------------------------------
      break;

    case BOTTOMRIGHT:
      psiKLeft  = psi[ sGrid.getKLeft ( i ) ];
      psiKAbove = psi[ sGrid.getKAbove ( i ) ];

      ALeft  = sGrid.getAxLeft ( i );
      AAbove = sGrid.getAyAbove ( i );

      res = ( - psiK      * 2.0
              + psiKLeft  * exp ( I*ALeft *h )
              + psiKAbove * exp ( -I*AAbove*h ) ) * I/ ( sqrt ( 2 ) *h );
      break;

    case TOPRIGHT:
      psiKLeft  = psi[ sGrid.getKLeft ( i ) ];
      psiKBelow = psi[ sGrid.getKBelow ( i ) ];

      ALeft  = sGrid.getAxLeft ( i );
      ABelow = sGrid.getAyBelow ( i );

      res = ( - psiK      * 2.0
              + psiKLeft  * exp ( I*ALeft *h )
              + psiKBelow * exp ( I*ABelow*h ) ) * I/ ( sqrt ( 2 ) *h );
      break;

    case TOPLEFT:
      psiKRight = psi[ sGrid.getKRight ( i ) ];
      psiKBelow = psi[ sGrid.getKBelow ( i ) ];

      ARight = sGrid.getAxRight ( i );
      ABelow = sGrid.getAyBelow ( i );

      res = ( - psiK      * 2.0
              + psiKRight * exp ( -I*ARight*h )
              + psiKBelow * exp ( I*ABelow*h ) ) * I/ ( sqrt ( 2 ) *h );
      break;

    case BOTTOM:
      // -------------------------------------------------------------------
      // normal derivative
      psiKAbove = psi[ sGrid.getKAbove ( i ) ];
      AAbove = sGrid.getAyAbove ( i );
      res = ( - psiK
              + psiKAbove * exp ( -I*AAbove*h ) ) * I/h;
      // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           // formulated with outer point, which is then eliminated by
//           // boundary condition
//           psiKLeft  = psi[ sGrid.getKLeft (i) ];
//           psiKRight = psi[ sGrid.getKRight(i) ];
//           psiKAbove = psi[ sGrid.getKAbove(i) ];
//
//           ALeft  = sGrid.getAxLeft ( i );
//           ARight = sGrid.getAxRight( i );
//           AAbove = sGrid.getAyAbove( i );
//           res = (   psiK*      (-3.0)
//                   + psiKLeft*  exp( I*ALeft *h) + psiKRight* exp(-I*ARight*h)
//                                                 + psiKAbove* exp(-I*AAbove*h) ) / (h*h)
//                 + psiK * (1-norm(psiK));
//           // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           // formulated with outer point, which is then eliminated by
//           // boundary condition
//           res = psiK;
//           // -------------------------------------------------------------------

      break;

    case RIGHT:
      // -------------------------------------------------------------------
      // normal derivative
      psiKLeft = psi[ sGrid.getKLeft ( i ) ];
      ALeft  = sGrid.getAxLeft ( i );
      res = ( - psiK
              + psiKLeft * exp ( I*ALeft*h ) ) * I/h;
      // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           // formulated with outer point, which is then eliminated by
//           // boundary condition
//           psiKLeft  = psi[ sGrid.getKLeft (i) ];
//           psiKBelow = psi[ sGrid.getKBelow(i) ];
//           psiKAbove = psi[ sGrid.getKAbove(i) ];
//
//           ALeft  = sGrid.getAxLeft( i );
//           ABelow = sGrid.getAyBelow( i );
//           AAbove = sGrid.getAyAbove( i );
//           res = (   psiK*      (-3.0)
//                   + psiKLeft*  exp( I*ALeft *h)
//                   + psiKBelow* exp( I*ABelow*h) + psiKAbove* exp(-I*AAbove*h) ) / (h*h)
//                 + psiK * (1-norm(psiK));
//           // -------------------------------------------------------------------
      break;

    case TOP:
      // -------------------------------------------------------------------
      // normal derivative
      psiKBelow = psi[ sGrid.getKBelow ( i ) ];
      ABelow = sGrid.getAyBelow ( i );
      res = ( - psiK
              + psiKBelow * exp ( I*ABelow*h ) ) * I/h;
      // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           // formulated with outer point, which is then eliminated by
//           // boundary condition
//           psiKLeft  = psi[ sGrid.getKLeft (i) ];
//           psiKRight = psi[ sGrid.getKRight(i) ];
//           psiKBelow = psi[ sGrid.getKBelow(i) ];
//
//           ALeft  = sGrid.getAxLeft ( i );
//           ARight = sGrid.getAxRight( i );
//           ABelow = sGrid.getAyBelow( i );
//           res = (   psiK*      (-3.0)
//                   + psiKLeft*  exp( I*ALeft *h) + psiKRight* exp(-I*ARight*h)
//                   + psiKBelow* exp( I*ABelow*h)                               ) / (h*h)
//                 + psiK * (1-norm(psiK));
//           // -------------------------------------------------------------------

      break;

    case LEFT:
      // -------------------------------------------------------------------
      // normal derivative
      psiKRight = psi[ sGrid.getKRight ( i ) ];
      ARight = sGrid.getAxRight ( i );
      res = ( - psiK
              + psiKRight * exp ( -I*ARight*h ) ) * I/h;
      // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           // formulated with outer point, which is then eliminated by
//           // boundary condition
//           psiKRight = psi[ sGrid.getKRight(i) ];
//           psiKBelow = psi[ sGrid.getKBelow(i) ];
//           psiKAbove = psi[ sGrid.getKAbove(i) ];
//
//           ARight = sGrid.getAxRight( i );
//           ABelow = sGrid.getAyBelow( i );
//           AAbove = sGrid.getAyAbove( i );
//           res = (   psiK*      (-3.0)
//                                                 + psiKRight* exp(-I*ARight*h)
//                   + psiKBelow* exp( I*ABelow*h) + psiKAbove* exp(-I*AAbove*h) ) / (h*h)
//                 + psiK * (1-norm(psiK));
//           // -------------------------------------------------------------------

      break;

    case INTERIOR:
      psiKLeft  = psi[ sGrid.getKLeft ( i ) ];
      psiKRight = psi[ sGrid.getKRight ( i ) ];
      psiKBelow = psi[ sGrid.getKBelow ( i ) ];
      psiKAbove = psi[ sGrid.getKAbove ( i ) ];

      ALeft  = sGrid.getAxLeft ( i );
      ARight = sGrid.getAxRight ( i );
      ABelow = sGrid.getAyBelow ( i );
      AAbove = sGrid.getAyAbove ( i );

      // -------------------------------------------------------------------
      res = ( psiK* ( -4.0 )
              + psiKLeft*  exp ( I*ALeft *h ) + psiKRight* exp ( -I*ARight*h )
              + psiKBelow* exp ( I*ABelow*h ) + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // -------------------------------------------------------------------
      break;

    default:
      throw glException ( "GinzburgLandau::computeGl",
                          "Illegal equationType "
                          + EpetraExt::toString ( eqType ) + "." );
    }

  // return the result
  return res;
}
// =============================================================================



// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
void GinzburgLandau::getJacobianRow ( const int                         eqnum,
                                      const std::vector<double_complex> &psi,
                                      std::vector<int>                  &columnIndicesPsi,
                                      std::vector<double_complex>       &valuesPsi,
                                      std::vector<int>                  &columnIndicesPsiConj,
                                      std::vector<double_complex>       &valuesPsiConj )
{
  computeJacobianRow ( VALUES,
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
void GinzburgLandau::getJacobianRowSparsity ( const int        eqnum,
    std::vector<int> &columnIndicesPsi,
    std::vector<int> &columnIndicesPsiConj )
{
  // create dummy arguments
  std::vector<double_complex> psi;
  std::vector<double_complex> valuesPsi,
  valuesPsiConj;

  computeJacobianRow ( SPARSITY,
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
void GinzburgLandau::computeJacobianRow ( const filltype                    ft,
    const int                         eqnum,
    const std::vector<double_complex> &psi,
    std::vector<int>                  &columnIndicesPsi,
    std::vector<double_complex>       &valuesPsi,
    std::vector<int>                  &columnIndicesPsiConj,
    std::vector<double_complex>       &valuesPsiConj )
{
  int          i[2];
  int k, kLeft, kRight, kBelow, kAbove;
  int numEntriesPsi, numEntriesPsiConj;
  double ARight, ALeft, AAbove, ABelow;
  double h = sGrid.getH();
  equationType eqType;

  // get the equation type from the running index
  getEquationType ( eqnum, eqType, i );

  // needed everywhere
  k = sGrid.i2k ( i );

  switch ( eqType )
    {
    case BOTTOMLEFT:
//           // -------------------------------------------------------------------
//           kRight = sGrid.getKRight( i );
//           kAbove = sGrid.getKAbove( i );
//
//           numEntriesPsi = 3;
//           columnIndicesPsi.resize(numEntriesPsi);
//           columnIndicesPsi[0] = k;
//           columnIndicesPsi[1] = kRight;
//           columnIndicesPsi[2] = kAbove;
//
//           if (ft==VALUES) {
//               ARight = sGrid.getAxRight( i );
//               AAbove = sGrid.getAyAbove( i );
//               valuesPsi.resize(numEntriesPsi);
//               valuesPsi[0] = -2.0             * I/(sqrt(2)*h);
//               valuesPsi[1] = exp(-I*ARight*h) * I/(sqrt(2)*h);
//               valuesPsi[2] = exp(-I*AAbove*h) * I/(sqrt(2)*h);
//           }
//           // -------------------------------------------------------------------


//       // -------------------------------------------------------------------
      kRight = sGrid.getKRight ( i );
//       kAbove = sGrid.getKAbove ( i );
// 
      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kRight;
//       columnIndicesPsi[2] = kAbove;

      if ( ft==VALUES )
        {
//           ARight = sGrid.getAxRight ( i );
//           AAbove = sGrid.getAyAbove ( i );

          valuesPsi.resize ( numEntriesPsi );
	  valuesPsi[0] = - 2.0 / (h*h)
	                 + ( 1 - 2.0*norm ( psi[k] ) );
          valuesPsi[1] = -30.0;
//           valuesPsi[2] = exp ( -I*AAbove*h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;
      if ( ft==VALUES )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psi[k]*psi[k];
        }

// std::cout << "--------------------------------------" << std::endl;
// for (unsigned int l=0; l<columnIndicesPsi.size(); l++ )
//     std::cout << " columnIndicesPsi[" << l << "] = " << columnIndicesPsi[l] << std::endl;
// std::cout << "--------------------------------------" << std::endl;
// for (unsigned int l=0; l<valuesPsi.size(); l++ )
//     std::cout << " valuesPsi[" << l << "] = " << valuesPsi[l] << std::endl;
// std::cout << "--------------------------------------" << std::endl;
// for (unsigned int l=0; l<columnIndicesPsiConj.size(); l++ )
//     std::cout << " columnIndicesPsiConj[" << l << "] = " << columnIndicesPsiConj[l] << std::endl;
// std::cout << "--------------------------------------" << std::endl;
// for (unsigned int l=0; l<valuesPsiConj.size(); l++ )
//     std::cout << " valuesPsiConj[" << l << "] = " << valuesPsiConj[l] << std::endl;
// std::cout << "--------------------------------------" << std::endl;

      // -------------------------------------------------------------------

      break;

    case BOTTOMRIGHT:
      kLeft  = sGrid.getKLeft ( i );
      kAbove = sGrid.getKAbove ( i );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kAbove;

      if ( ft==VALUES )
        {
          ALeft    = sGrid.getAxLeft ( i );
          AAbove   = sGrid.getAyAbove ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -2.0                * I/ ( sqrt ( 2 ) *h );
          valuesPsi[1] = exp (  I*ALeft *h ) * I/ ( sqrt ( 2 ) *h );
          valuesPsi[2] = exp ( -I*AAbove*h ) * I/ ( sqrt ( 2 ) *h );
        }

      break;

    case TOPRIGHT:
      kLeft  = sGrid.getKLeft ( i );
      kBelow = sGrid.getKBelow ( i );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kBelow;

      if ( ft==VALUES )
        {
          ALeft    = sGrid.getAxLeft ( i );
          ABelow   = sGrid.getAyBelow ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -2.0             * I/ ( sqrt ( 2 ) *h );
          valuesPsi[1] = exp ( I*ALeft *h ) * I/ ( sqrt ( 2 ) *h );
          valuesPsi[2] = exp ( I*ABelow*h ) * I/ ( sqrt ( 2 ) *h );
        }

      break;

    case TOPLEFT:
      kRight = sGrid.getKRight ( i );
      kBelow = sGrid.getKBelow ( i );

      numEntriesPsi = 3;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kRight;
      columnIndicesPsi[2] = kBelow;

      if ( ft==VALUES )
        {
          ARight    = sGrid.getAxRight ( i );
          ABelow    = sGrid.getAyBelow ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -2.0             * I/ ( sqrt ( 2 ) *h );
          valuesPsi[1] = exp ( -I*ARight*h ) * I/ ( sqrt ( 2 ) *h );
          valuesPsi[2] = exp ( I*ABelow*h ) * I/ ( sqrt ( 2 ) *h );
        }

      break;

    case BOTTOM:
      // -------------------------------------------------------------------
      // normal derivative
      kAbove = sGrid.getKAbove ( i );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kAbove;

      if ( ft==VALUES )
        {
          AAbove = sGrid.getAyAbove ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -1.0             * I/h;
          valuesPsi[1] = exp ( -I*AAbove*h ) * I/h;
        }
      // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           // formulated with outer point, which is then eliminated by
//           // boundary condition
//           kLeft  = sGrid.getKLeft ( i );
//           kRight = sGrid.getKRight( i );
//           kAbove = sGrid.getKAbove( i );
//
//           numEntriesPsi = 4;
//           columnIndicesPsi.resize(numEntriesPsi);
//           columnIndicesPsi[0] = k;
//           columnIndicesPsi[1] = kLeft;
//           columnIndicesPsi[2] = kRight;
//           columnIndicesPsi[3] = kAbove;
//
//           if (ft==VALUES) {
//               ALeft  = sGrid.getAxLeft ( i );
//               ARight = sGrid.getAxRight( i );
//               AAbove = sGrid.getAyAbove( i );
//               valuesPsi.resize(numEntriesPsi);
//               valuesPsi[0] = - 3.0             /(h*h)
//                              + (1 - 2.0*norm(psi[k]));
//               valuesPsi[1] = exp( I*ALeft *h) /(h*h);
//               valuesPsi[2] = exp(-I*ARight*h) /(h*h);
//               valuesPsi[3] = exp(-I*AAbove*h) /(h*h);
//           }
//
//           numEntriesPsiConj = 1;
//           columnIndicesPsiConj.resize(numEntriesPsiConj);
//           columnIndicesPsiConj[0] = k;
//
//           if (ft==VALUES) {
//               valuesPsiConj.resize(numEntriesPsiConj);
//               valuesPsiConj[0] = -psi[k]*psi[k];
//           }
//           // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           numEntriesPsiConj = 1;
//           columnIndicesPsiConj[0] = k;
//
//           if (ft==VALUES) {
//               valuesPsi.resize(numEntriesPsiConj);
//               valuesPsi[0] = 1.0;
//           }
//           // -------------------------------------------------------------------
      break;

    case RIGHT:
      // -------------------------------------------------------------------
      // normal derivative
      kLeft = sGrid.getKLeft ( i );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;

      if ( ft==VALUES )
        {
          ALeft = sGrid.getAxLeft ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -1.0            * I/h;
          valuesPsi[1] = exp ( I*ALeft*h ) * I/h;
        }
      // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           // formulated with outer point, which is then eliminated by
//           // boundary condition
//           kBelow = sGrid.getKBelow( i );
//           kAbove = sGrid.getKAbove( i );
//           kLeft  = sGrid.getKLeft ( i );
//
//           numEntriesPsi = 4;
//           columnIndicesPsi.resize(numEntriesPsi);
//           columnIndicesPsi[0] = k;
//           columnIndicesPsi[1] = kBelow;
//           columnIndicesPsi[2] = kAbove;
//           columnIndicesPsi[3] = kLeft;
//
//           if (ft==VALUES) {
//               ABelow = sGrid.getAyBelow( i );
//               AAbove = sGrid.getAyAbove( i );
//               ALeft  = sGrid.getAxLeft ( i );
//               valuesPsi.resize(numEntriesPsi);
//               valuesPsi[0] = -3.0             /(h*h);
//               valuesPsi[1] = exp( I*ABelow*h) /(h*h);
//               valuesPsi[2] = exp(-I*AAbove*h) /(h*h);
//               valuesPsi[3] = exp( I*ALeft *h) /(h*h);
//           }
//
//           numEntriesPsiConj = 1;
//           columnIndicesPsiConj.resize(numEntriesPsiConj);
//           columnIndicesPsiConj[0] = k;
//
//           if (ft==VALUES) {
//               valuesPsiConj.resize(numEntriesPsiConj);
//               valuesPsiConj[0] = -psi[k]*psi[k];
//           }
//           // -------------------------------------------------------------------
      break;

    case TOP:
      // -------------------------------------------------------------------
      // normal derivative
      kBelow = sGrid.getKBelow ( i );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kBelow;

      if ( ft==VALUES )
        {
          ABelow = sGrid.getAyBelow ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -1.0             * I/h;
          valuesPsi[1] = exp ( I*ABelow*h ) * I/h;
        }
      // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           // formulated with outer point, which is then eliminated by
//           // boundary condition
//           kBelow = sGrid.getKBelow( i );
//           kRight = sGrid.getKRight( i );
//           kLeft  = sGrid.getKLeft ( i );
//
//           numEntriesPsi = 4;
//           columnIndicesPsi.resize(numEntriesPsi);
//           columnIndicesPsi[0] = k;
//           columnIndicesPsi[1] = kBelow;
//           columnIndicesPsi[2] = kLeft;
//           columnIndicesPsi[3] = kRight;
//
//           if (ft==VALUES) {
//               ABelow = sGrid.getAyBelow( i );
//               ALeft  = sGrid.getAxLeft ( i );
//               ARight = sGrid.getAxRight( i );
//               valuesPsi.resize(numEntriesPsi);
//               valuesPsi[0] = -3.0             /(h*h);
//               valuesPsi[1] = exp( I*ABelow*h) /(h*h);
//               valuesPsi[2] = exp( I*ALeft *h) /(h*h);
//               valuesPsi[3] = exp(-I*ARight*h) /(h*h);
//           }
//
//           numEntriesPsiConj = 1;
//           columnIndicesPsiConj.resize(numEntriesPsiConj);
//           columnIndicesPsiConj[0] = k;
//
//           if (ft==VALUES) {
//               valuesPsiConj.resize(numEntriesPsiConj);
//               valuesPsiConj[0] = -psi[k]*psi[k];
//           }
//           // -------------------------------------------------------------------

      break;

    case LEFT:
      // -------------------------------------------------------------------
      // normal derivative
      kRight = sGrid.getKRight ( i );

      numEntriesPsi = 2;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kRight;

      if ( ft==VALUES )
        {
          ARight = sGrid.getAxRight ( i );
          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = -1.0             * I/h;
          valuesPsi[1] = exp ( -I*ARight*h ) * I/h;
        }
      // -------------------------------------------------------------------

//           // -------------------------------------------------------------------
//           // formulated with outer point, which is then eliminated by
//           // boundary condition
//           kBelow = sGrid.getKBelow( i );
//           kAbove = sGrid.getKAbove( i );
//           kRight = sGrid.getKRight( i );
//
//           numEntriesPsi = 4;
//           columnIndicesPsi.resize(numEntriesPsi);
//           columnIndicesPsi[0] = k;
//           columnIndicesPsi[1] = kBelow;
//           columnIndicesPsi[2] = kAbove;
//           columnIndicesPsi[3] = kRight;
//
//           if (ft==VALUES) {
//               ABelow = sGrid.getAyBelow( i );
//               AAbove = sGrid.getAyAbove( i );
//               ARight = sGrid.getAxRight( i );
//               valuesPsi.resize(numEntriesPsi);
//               valuesPsi[0] = -3.0             /(h*h);
//               valuesPsi[1] = exp( I*ABelow*h) /(h*h);
//               valuesPsi[2] = exp(-I*AAbove*h) /(h*h);
//               valuesPsi[3] = exp(-I*ARight*h) /(h*h);
//           }
//
//           numEntriesPsiConj = 1;
//           columnIndicesPsiConj.resize(numEntriesPsiConj);
//           columnIndicesPsiConj[0] = k;
//
//           if (ft==VALUES) {
//               valuesPsiConj.resize(numEntriesPsiConj);
//               valuesPsiConj[0] = -psi[k]*psi[k];
//           }
//           // -------------------------------------------------------------------
      break;

    case INTERIOR:
      kRight = sGrid.getKRight ( i );
      kLeft  = sGrid.getKLeft ( i );
      kAbove = sGrid.getKAbove ( i );
      kBelow = sGrid.getKBelow ( i );

      numEntriesPsi = 5;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kRight;
      columnIndicesPsi[3] = kBelow;
      columnIndicesPsi[4] = kAbove;

      if ( ft==VALUES )
        {
          ALeft  = sGrid.getAxLeft ( i );
          ARight = sGrid.getAxRight ( i );
          ABelow = sGrid.getAyBelow ( i );
          AAbove = sGrid.getAyAbove ( i );

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = - 4.0            / ( h*h )
                         + ( 1 - 2.0*norm ( psi[k] ) );
          valuesPsi[1] = exp (  I*ALeft *h ) / ( h*h );
          valuesPsi[2] = exp ( -I*ARight*h ) / ( h*h );
          valuesPsi[3] = exp (  I*ABelow*h ) / ( h*h );
          valuesPsi[4] = exp ( -I*AAbove*h ) / ( h*h );
        }

      numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;
      if ( ft==VALUES )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psi[k]*psi[k];
        }
      break;

    default:
      throw glException ( "GinzburgLandau::getJacobianRow",
                          "Illegal equationType"
                          + EpetraExt::toString ( eqType ) + "." );
    }
}
// =============================================================================


// =============================================================================
// calculate the free energy of a state
double GinzburgLandau::freeEnergy ( const std::vector<double_complex> &psi )
{
  double                  energy = 0.0,
                                   h      = sGrid.getH();
  StaggeredGrid::nodeType nt;

  for ( unsigned int k=0; k<psi.size(); k++ )
    {
      nt = sGrid.k2nodeType ( k );
      if ( nt==StaggeredGrid::CORNER )
        energy -=       h*h * pow ( norm ( psi[k] ),2 );
      else if ( nt==StaggeredGrid::EDGE )
        energy -= 0.5*  h*h * pow ( norm ( psi[k] ),2 );
      else if ( nt==StaggeredGrid::INTERIOR )
        energy -= 0.25* h*h * pow ( norm ( psi[k] ),2 );
      else
        {
          std::string message = "Illegal node type "
                                + EpetraExt::toString ( nt )
                                + ".";
          throw glException ( "GinzburgLandau::freeEnergy",
                              message );
        }
    }

  return energy;
}
// =============================================================================
