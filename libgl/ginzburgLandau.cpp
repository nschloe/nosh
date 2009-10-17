#include "ginzburgLandau.h"

#include "glBoundaryConditionsVirtual.h"
#include "glException.h"

#include <iostream>

#include <fstream> // for the VTK writer
#include <string>

#include <vector>

// #include <EpetraExt_HDF5.h> // for output in XDMF format
// #include <Epetra_MultiVector.h>
// #include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>


// really needed?
// --> reduceAllAndScatter in freeEnergy()
#include <Teuchos_Comm.hpp>

// #include <Tpetra_Vector.hpp>

#include <EpetraExt_Utils.h> // for the toString function

// complex unit
const double_complex I ( 0,1 );

// =============================================================================
// Class constructor
GinzburgLandau::GinzburgLandau ( int    nx,
                                 double edgelength,
                                 double h0,
                                 Teuchos::RCP<GlBoundaryConditionsVirtual> bc
                               ) :
    sGrid ( StaggeredGrid::StaggeredGrid ( nx,
                                           edgelength,
                                           h0 ) ),
    boundaryConditions_( bc )
{
}
// =============================================================================
// Destructor
GinzburgLandau::~GinzburgLandau()
{
}
// =============================================================================
StaggeredGrid::StaggeredGrid* GinzburgLandau::getStaggeredGrid()
{
  return &sGrid;
}
// =============================================================================
// Defines a mapping of all the GL equations to a running index.
// The equations in the GL context are each associated with and
// uniquely identified by a particular node
// (around which the equation is centered, e.g., the center node with the five-
// point stencil).
void GinzburgLandau::getEquationType ( const int           eqnum,
                                       equationType        &eqType,
                                       int                 &eqIndex )
{
  int Nx = sGrid.getNx();
  int numBoundaryEquations = 4*Nx;
  int numTotalEquations =  (Nx+1) * (Nx+1);
  
  if ( eqnum<numBoundaryEquations )
    {
      eqType = BOUNDARY;
      eqIndex = eqnum;
    }
  else if ( eqnum<numTotalEquations )
    {
      eqType = INTERIOR;
      eqIndex = eqnum - numBoundaryEquations;
    }
  else
    {
      throw glException ( "GinzburgLandau::getEquationType",
                          "Illegal running index eqnum="
                          + EpetraExt::toString ( eqnum ) );
    }
}
// =============================================================================
Tpetra::MultiVector<double_complex,int>
GinzburgLandau::computeGlVector( Tpetra::MultiVector<double_complex,int> psi )
{
  // setup output vector with the same map as psi
  Tpetra::MultiVector<double_complex,int> glVec( psi.getMap(), 1, true );
  
  for ( unsigned int k=0; k<psi.getLocalLength(); k++ ) {
    int            globalIndex = psi.getMap()->getGlobalElement(k);
    double_complex z           = computeGl( globalIndex, psi );
    glVec.replaceLocalValue( k, 0, z );
  }
  
  return glVec;
}
// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
double_complex
GinzburgLandau::computeGl ( const int                                     eqnum,
                            const Tpetra::MultiVector<double_complex,int> &psi   )
{
  // the preliminary result type
  double_complex res;
  equationType eqType;
  int eqIndex;

  // get the equation type and the subindex from the running index
  getEquationType ( eqnum, eqType, eqIndex );

  switch ( eqType )
    {
    // -------------------------------------------------------------------------
    case BOUNDARY:
      res = boundaryConditions_->getGlEntry( eqIndex, psi, sGrid );
      break;
    // -------------------------------------------------------------------------
    case INTERIOR: {
      // translate eqIndex to i
      // TODO: Unify this with what happens at the Jacobian.
      int Nx = sGrid.getNx();
      Teuchos::Array<int> i(2);
      i[0] = eqIndex % ( Nx-1 ) + 1;
      i[1] = eqIndex / ( Nx-1 ) + 1;
  
      Teuchos::ArrayRCP<const double_complex> psiView =
                                                  psi.getVector(0)->get1dView();

      double_complex psiK      = psiView[ sGrid.i2k       ( i ) ];
      double_complex psiKLeft  = psiView[ sGrid.getKLeft  ( i ) ];
      double_complex psiKRight = psiView[ sGrid.getKRight ( i ) ];
      double_complex psiKBelow = psiView[ sGrid.getKBelow ( i ) ];
      double_complex psiKAbove = psiView[ sGrid.getKAbove ( i ) ];

      double ALeft  = sGrid.getAxLeft  ( i );
      double ARight = sGrid.getAxRight ( i );
      double ABelow = sGrid.getAyBelow ( i );
      double AAbove = sGrid.getAyAbove ( i );

      double h = sGrid.getH();

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      res = ( psiK* ( -4.0 )
              + psiKLeft*  exp ( I*ALeft *h ) + psiKRight* exp ( -I*ARight*h )
              + psiKBelow* exp ( I*ABelow*h ) + psiKAbove* exp ( -I*AAbove*h ) ) / ( h*h )
            + psiK * ( 1-norm ( psiK ) );
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
      break;
    // -------------------------------------------------------------------------
    default:
      throw glException ( "GinzburgLandau::computeGl",
                          "Illegal equationType "
                          + EpetraExt::toString ( eqType ) + "." );
    }

  // return the result
  return res;
}
// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
void GinzburgLandau::getJacobianRow ( const int                                     eqnum,
                                      const Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi,
                                      std::vector<int>                              &columnIndicesPsi,
                                      std::vector<double_complex>                   &valuesPsi,
                                      std::vector<int>                              &columnIndicesPsiConj,
                                      std::vector<double_complex>                   &valuesPsiConj )
{
  computeJacobianRow ( true,
                       eqnum,
                       psi,
                       columnIndicesPsi,
                       valuesPsi,
                       columnIndicesPsiConj,
                       valuesPsiConj );

  return;
}
// =============================================================================
void GinzburgLandau::getJacobianRowSparsity ( const int        eqnum,
                                              std::vector<int> &columnIndicesPsi,
                                              std::vector<int> &columnIndicesPsiConj )
{
  // create dummy arguments
  Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi = Teuchos::ENull();

  std::vector<double_complex> valuesPsi, valuesPsiConj;

  computeJacobianRow ( false,
                       eqnum,
                       psi,
                       columnIndicesPsi,
                       valuesPsi,
                       columnIndicesPsiConj,
                       valuesPsiConj );

  return;
}
// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
void GinzburgLandau::computeJacobianRow ( const bool                                    fillValues,
                                          const int                                     eqnum,
                                          const Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi,
                                          std::vector<int>                              &columnIndicesPsi,
                                          std::vector<double_complex>                   &valuesPsi,
                                          std::vector<int>                              &columnIndicesPsiConj,
                                          std::vector<double_complex>                   &valuesPsiConj )
{ 
  equationType eqType;
  int          eqIndex;

  // get the equation type and the subindex from the running index
  getEquationType ( eqnum, eqType, eqIndex );
  
  switch ( eqType )
    {
    case BOUNDARY:
        // ---------------------------------------------------------------------
        boundaryConditions_->getGlJacobianRow ( eqIndex,
                                                psi,
                                                sGrid,
                                                fillValues,
                                                columnIndicesPsi,
                                                valuesPsi,
                                                columnIndicesPsiConj,
                                                valuesPsiConj );
        // ---------------------------------------------------------------------
        break;

    case INTERIOR: {
      // ---------------------------------------------------------------------
      // translate eqIndex to i
      // TODO: Unify this with what happens at the Jacobian.
      int Nx = sGrid.getNx();
      Teuchos::Array<int> i(2);
      i[0] = eqIndex % ( Nx-1 ) + 1;
      i[1] = eqIndex / ( Nx-1 ) + 1;

      int k      = sGrid.i2k( i );
      int kRight = sGrid.getKRight ( i );
      int kLeft  = sGrid.getKLeft ( i );
      int kAbove = sGrid.getKAbove ( i );
      int kBelow = sGrid.getKBelow ( i );

      int numEntriesPsi = 5;
      columnIndicesPsi.resize ( numEntriesPsi );
      columnIndicesPsi[0] = k;
      columnIndicesPsi[1] = kLeft;
      columnIndicesPsi[2] = kRight;
      columnIndicesPsi[3] = kBelow;
      columnIndicesPsi[4] = kAbove;

      Teuchos::ArrayRCP<const double_complex> psiView;
      if ( fillValues )
        {
	  psiView = psi->getVector(0)->get1dView();

          double ALeft  = sGrid.getAxLeft ( i );
          double ARight = sGrid.getAxRight ( i );
          double ABelow = sGrid.getAyBelow ( i );
          double AAbove = sGrid.getAyAbove ( i );

          double h = sGrid.getH();

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = - 4.0            / ( h*h )
                         + ( 1 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp (  I*ALeft *h ) / ( h*h );
          valuesPsi[2] = exp ( -I*ARight*h ) / ( h*h );
          valuesPsi[3] = exp (  I*ABelow*h ) / ( h*h );
          valuesPsi[4] = exp ( -I*AAbove*h ) / ( h*h );
        }

      int numEntriesPsiConj = 1;
      columnIndicesPsiConj.resize ( numEntriesPsiConj );
      columnIndicesPsiConj[0] = k;
      if ( fillValues )
        {
          valuesPsiConj.resize ( numEntriesPsiConj );
          valuesPsiConj[0] = -psiView[k]*psiView[k];
        }
    }
      break;

    default:
      throw glException ( "GinzburgLandau::getJacobianRow",
                          "Illegal equationType"
                          + EpetraExt::toString ( eqType ) + "." );
    }
}
// =============================================================================
// calculate the free energy of a state
double GinzburgLandau::freeEnergy ( const Tpetra::MultiVector<double_complex,int> &psi )
{
  double localEnergy = 0.0;
  double h      = sGrid.getH();
  StaggeredGrid::nodeType nt;
  
  Teuchos::ArrayRCP<const double_complex> psiView = psi.getVector(0)->get1dView();
  
  // sum up the energy on each processor
  for ( int k=0; k<psi.getLocalLength(); k++ )
    {
      nt = sGrid.k2nodeType ( psi.getMap()->getGlobalElement(k) );
      if ( nt==StaggeredGrid::CORNER )
        localEnergy -= 0.25* h*h * pow ( norm ( psiView[k] ),2 );
      else if ( nt==StaggeredGrid::EDGE )
        localEnergy -= 0.5*  h*h * pow ( norm ( psiView[k] ),2 );
      else if ( nt==StaggeredGrid::INTERIOR )
        localEnergy -=       h*h * pow ( norm ( psiView[k] ),2 );
      else
        {
          std::string message = "Illegal node type "
                                + EpetraExt::toString ( nt )
                                + ".";
          throw glException ( "GinzburgLandau::freeEnergy",
                              message );
        }
    }

    // reduce and scatter such that energy is available on
    // all cores
    int count = 1; // send *one* integer
    int numProcs =  psi.getMap()->getComm()->getSize();
    Teuchos::Array<double> sendBuff(count), recvBuff(count);
    
    // fill send buffer
    sendBuff[0] = localEnergy;
    
    Teuchos::Array<int> recvCounts(numProcs);
// fill recvCounts with {1,...,1}
  int numItemsPerProcess = 1;
  std::fill(recvCounts.begin(), recvCounts.end(), numItemsPerProcess);
  
    Teuchos::reduceAllAndScatter( *(psi.getMap()->getComm()),
                                  Teuchos::REDUCE_SUM,
                                  count,
                                  &sendBuff[0],
                                  &recvCounts[0],
                                  &recvBuff[0]
                                );

    int globalEnergy = recvBuff[0];

  return globalEnergy;
}
// =============================================================================
// count the number of vortices by the total phase change along the boundary
// TODO:
// make this work in multicore environments
int GinzburgLandau::countVortices ( const Tpetra::MultiVector<double_complex,int> &psi )
{ 
  // this function only works 
  int numProcs = psi.getMap()->getComm()->getSize();
  if ( numProcs!=1 )
    return -1;
  
  int numVortices = 0;
  Teuchos::Array<int> i(2,0.0);
  int k;
  int Nx = sGrid.getNx();
  
  const double pi = 3.14159265358979323846264338327950288419716939937510;
  const double threshold = 1.5*pi; // Consider jumps in the argument greater
                                   // than this phase jumps.

  // Get a view of the whole vector.
  // Remember: This only works with one core.
  Teuchos::ArrayRCP<const double_complex> psiView = psi.getVector(0)->get1dView();

    cout << "ccc" << endl;
  
  // origin -- our first index
  k = sGrid.i2k( i );

  double angle = arg(psiView[k]);
  double anglePrev;

  // lower border
  i[1] = 0;
  for ( int l=1; l<Nx+1; l++ ) {
      anglePrev = angle;
      i[0] = l;
      k = sGrid.i2k( i );
      angle = arg(psiView[k]);
      if ( abs(angle-anglePrev)>threshold )
          numVortices++;
  }

  // right border
  i[0] = Nx;
  for ( int l=1; l<Nx+1; l++ ) {
      anglePrev = angle;
      i[1] = l;
      k = sGrid.i2k( i );
      angle = arg(psiView[k]);
      if ( abs(angle-anglePrev)>threshold )
          numVortices++;
  }

  // top border
  i[1] = Nx;
  for ( int l=1; l<Nx+1; l++ ) {
      anglePrev = angle;
      i[0] = Nx-l;
      k = sGrid.i2k( i );
      angle = arg(psiView[k]);
      if ( abs(angle-anglePrev)>threshold )
          numVortices++;
  }

  // left border
  i[0] = 0;
  for ( int l=1; l<Nx+1; l++ ) {
      anglePrev = angle;
      i[1] = Nx-l;
      k = sGrid.i2k( i );
      angle = arg(psiView[k]);
      if ( abs(angle-anglePrev)>threshold )
          numVortices++;
  }

  return numVortices;
}
// =============================================================================