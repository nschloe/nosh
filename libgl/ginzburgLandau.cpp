#include "ginzburgLandau.h"

#include "glBoundaryConditionsVirtual.h"
#include "glException.h"

#include <Teuchos_RCP.hpp>

// really needed?
// --> reduceAllAndScatter in freeEnergy()
#include <Teuchos_Comm.hpp>

#include <EpetraExt_Utils.h> // for the toString function

// complex unit
const double_complex I ( 0,1 );

// =============================================================================
// Class constructor
GinzburgLandau::GinzburgLandau ( Teuchos::RCP<Grid>                        &grid,
                                 Teuchos::RCP<MagneticVectorPotential>     &A,
                                 Teuchos::RCP<GlBoundaryConditionsVirtual> &bc
                               ) :
    grid_ ( grid ),
    A_ ( A ),
    boundaryConditions_ ( bc )
{
}
// =============================================================================
// Destructor
GinzburgLandau::~GinzburgLandau()
{
}
// =============================================================================
Teuchos::RCP<MagneticVectorPotential>
GinzburgLandau::getMagneticVectorPotential() const
{
  return A_;
}
// =============================================================================
Teuchos::RCP<Grid>
GinzburgLandau::getGrid() const
{
  return grid_;
}
// =============================================================================
// Defines a mapping of all the GL equations to a running index.
// The equations in the GL context are each associated with and
// uniquely identified by a particular node
// (around which the equation is centered, e.g., the center node with the five-
// point stencil).
void GinzburgLandau::getEquationType ( const int           eqnum,
                                       equationType        &eqType,
                                       int                 &eqIndex
                                     ) const
{
  int Nx = grid_->getNx();
  int numBoundaryEquations = 4*Nx;
  int numTotalEquations = ( Nx+1 ) * ( Nx+1 );

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
Tpetra::Vector<double_complex,int>
GinzburgLandau::computeGlVector ( const ComplexVector psi
                                ) const
{
  // setup output vector with the same map as psi
  ComplexVector glVec ( psi.getMap(), true );

  for ( unsigned int k=0; k<psi.getLocalLength(); k++ )
    {
      int            globalIndex = psi.getMap()->getGlobalElement ( k );
      double_complex z           = computeGl ( globalIndex, psi );
      glVec.replaceLocalValue ( k, z );
    }

  return glVec;
}
// =============================================================================
// Evaluate GL at the boundary node.
// Return value for equation #k.
double_complex
GinzburgLandau::computeGl ( const int           eqnum,
                            const ComplexVector &psi
                          ) const
{
  // the preliminary result type
  double_complex res;
  equationType eqType;
  int eqIndex;

  // get the equation type and the sub-index from the running index (0,1,2,3,...
  // of the boundary; 0,1,2,3,... of the interior)
  getEquationType ( eqnum, eqType, eqIndex );

  switch ( eqType )
    {
      // -------------------------------------------------------------------------
    case BOUNDARY:
      res = boundaryConditions_->getGlEntry ( eqIndex, psi, *grid_, *A_ );
      break;
      // -------------------------------------------------------------------------
    case INTERIOR:
    {
      // Translate eqIndex to k.
      // TODO: Unify this with what happens at the Jacobian.
      // The equations here are "centered" around a node. Take node 'k'
      // for equation 'eqIndex'.
      int k = eqnum;

      Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

      double_complex psiK      = psiView[ k ];
      double_complex psiKLeft  = psiView[ grid_->getKLeft ( k ) ];
      double_complex psiKRight = psiView[ grid_->getKRight( k ) ];
      double_complex psiKBelow = psiView[ grid_->getKBelow( k ) ];
      double_complex psiKAbove = psiView[ grid_->getKAbove( k ) ];

      Teuchos::RCP<Teuchos::Array<double> > xLeft ( grid_->getXLeft ( k ) );
      Teuchos::RCP<Teuchos::Array<double> > xRight( grid_->getXRight( k ) );
      Teuchos::RCP<Teuchos::Array<double> > xBelow( grid_->getXBelow( k ) );
      Teuchos::RCP<Teuchos::Array<double> > xAbove( grid_->getXAbove( k ) );

      double ALeft  = A_->getAx ( *xLeft  );
      double ARight = A_->getAx ( *xRight );
      double ABelow = A_->getAy ( *xBelow );
      double AAbove = A_->getAy ( *xAbove );

      double h = grid_->getH();

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
void GinzburgLandau::getJacobianRow ( const int                         eqnum,
                                      const Teuchos::RCP<ComplexVector> &psi,
                                      std::vector<int>                  &columnIndicesPsi,
                                      std::vector<double_complex>       &valuesPsi,
                                      std::vector<int>                  &columnIndicesPsiConj,
                                      std::vector<double_complex>       &valuesPsiConj
                                    ) const
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
                                              std::vector<int> &columnIndicesPsiConj
                                            ) const
{
  // create dummy arguments
  Teuchos::RCP<ComplexVector > psi = Teuchos::ENull();

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
void GinzburgLandau::computeJacobianRow ( const bool                        fillValues,
                                          const int                         eqnum,
                                          const Teuchos::RCP<ComplexVector> &psi,
                                          std::vector<int>                  &columnIndicesPsi,
                                          std::vector<double_complex>       &valuesPsi,
                                          std::vector<int>                  &columnIndicesPsiConj,
                                          std::vector<double_complex>       &valuesPsiConj
                                        ) const
{
  equationType eqType;
  int          eqIndex;

  // get the equation type and the sub-index from the running index
  getEquationType ( eqnum, eqType, eqIndex );

  switch ( eqType )
    {
    case BOUNDARY:
      // ---------------------------------------------------------------------
      boundaryConditions_->getGlJacobianRow ( eqIndex,
                                              psi,
                                              *grid_,
                                              *A_,
                                              fillValues,
                                              columnIndicesPsi,
                                              valuesPsi,
                                              columnIndicesPsiConj,
                                              valuesPsiConj );
      // ---------------------------------------------------------------------
      break;

    case INTERIOR:
    {
      // ---------------------------------------------------------------------
      // Translate eqIndex to k.
      // TODO: Unify this with what happens at the Jacobian.
      // The equations here are "centered" around a node. Take node 'k'
      // for equation 'eqIndex'.
      int k = eqnum;

      int kRight = grid_->getKRight ( k );
      int kLeft  = grid_->getKLeft  ( k );
      int kAbove = grid_->getKAbove ( k );
      int kBelow = grid_->getKBelow ( k );

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
          psiView = psi->get1dView();

          Teuchos::RCP<Teuchos::Array<double> > xLeft ( grid_->getXLeft ( k ) );
          Teuchos::RCP<Teuchos::Array<double> > xRight( grid_->getXRight( k ) );
          Teuchos::RCP<Teuchos::Array<double> > xBelow( grid_->getXBelow( k ) );
          Teuchos::RCP<Teuchos::Array<double> > xAbove( grid_->getXAbove( k ) );
          double ALeft  = A_->getAx ( *xLeft  );
          double ARight = A_->getAx ( *xRight );
          double ABelow = A_->getAy ( *xBelow );
          double AAbove = A_->getAy ( *xAbove );

          double h = grid_->getH();

          valuesPsi.resize ( numEntriesPsi );
          valuesPsi[0] = - 4.0            / ( h*h )
                         + ( 1 - 2.0*norm ( psiView[k] ) );
          valuesPsi[1] = exp ( I*ALeft *h ) / ( h*h );
          valuesPsi[2] = exp ( -I*ARight*h ) / ( h*h );
          valuesPsi[3] = exp ( I*ABelow*h ) / ( h*h );
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
double GinzburgLandau::freeEnergy ( const ComplexVector &psi
                                  ) const
{
  double localEnergy = 0.0;

  Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

  // sum up the energy on each processor
  for ( int k=0; k<psi.getLocalLength(); k++ )
    {
	  int kGlobal = psi.getMap()->getGlobalElement ( k );
	  double area = grid_->cellArea( kGlobal );
	  localEnergy -= area * pow ( norm ( psiView[k] ),2 );
    }
  
  // reduce and scatter such that energy is available on
  // all cores
  int count = 1; // send *one* integer
  int numProcs =  psi.getMap()->getComm()->getSize();
  Teuchos::Array<double> sendBuff ( count ), recvBuff ( count );

  // fill send buffer
  sendBuff[0] = localEnergy;
  
  Teuchos::Array<int> recvCounts ( numProcs );
// fill recvCounts with {1,...,1}
  int numItemsPerProcess = 1;
  std::fill ( recvCounts.begin(), recvCounts.end(), numItemsPerProcess );

  Teuchos::reduceAllAndScatter ( * ( psi.getMap()->getComm() ),
                                 Teuchos::REDUCE_SUM,
                                 count,
                                 &sendBuff[0],
                                 &recvCounts[0],
                                 &recvBuff[0]
                               );

  // normalize
  double maxEnergy = 1.0 * grid_->getGridDomainArea();
  double globalEnergy = recvBuff[0] / maxEnergy;

  return globalEnergy;
}
// =============================================================================
// Count the number of vortices by the total phase change along the boundary
// of the domain.
// TODO:
// Make this work in multicore environments.
int GinzburgLandau::getVorticity ( const ComplexVector &psi
                                 ) const
{
  // this function only works
  int numProcs = psi.getMap()->getComm()->getSize();
  if ( numProcs!=1 )
    return -1;

  int vorticity = 0;

  const double PI = 3.14159265358979323846264338327950288419716939937510;
  const double threshold = 1.5*PI; // Consider jumps in the argument greater
  // than this phase jumps.

  // Get a view of the whole vector.
  // Remember: This only works with one core.
  Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

  int numBorderPoints = grid_->getNumBorderPoints();

  int l = 0;
  double angle = arg( psiView[ grid_->borderNode(l) ] );
  double anglePrevious;
  for ( l=1; l<numBorderPoints; l++ ){
	  anglePrevious = angle;
	  angle = arg( psiView[grid_->borderNode(l)] );
      if ( angle-anglePrevious<-threshold )
        vorticity++;
      else if ( angle-anglePrevious>threshold )
        vorticity--;
  }

  // close the circle
  l = 0;
  anglePrevious = angle;
  angle = arg( psiView[grid_->borderNode(l)] );
  if ( angle-anglePrevious<-threshold )
    vorticity++;
  else if ( angle-anglePrevious>threshold )
    vorticity--;

  return vorticity;
}
// =============================================================================
