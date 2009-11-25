#include "ginzburgLandau.h"

#include "glBoundaryConditionsVirtual.h"

#include <Teuchos_RCP.hpp>

// really needed?
// --> reduceAllAndScatter in freeEnergy()
#include <Teuchos_Comm.hpp>

#include <Tpetra_Map.hpp>

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
void
GinzburgLandau::setH0(const double h0)
{
  A_->setH0( h0 );
}
// =============================================================================
void
GinzburgLandau::setEdgeLength( const double edgeLength)
{
  grid_->setEdgeLength( edgeLength );
}
// =============================================================================
int
GinzburgLandau::getNumUnknowns() const
{
  return grid_->getNumGridPoints();
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
  int numBoundaryEquations = grid_->getNumBoundaryPoints();
  int numTotalEquations = grid_->getNumGridPoints();

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
      TEST_FOR_EXCEPTION( true,
  			              std::logic_error,
  			              "Illegal running index   eqnum = " <<  eqnum );
    }
}
// =============================================================================
ComplexVector
GinzburgLandau::computeGlVector ( const Teuchos::RCP<ComplexVector> & psi
                                ) const
{
  // setup output vector with the same map as psi
  ComplexVector glVec ( psi->getMap(), true );

  for ( unsigned int k=0; k<psi->getLocalLength(); k++ )
    {
      int            globalIndex = psi->getMap()->getGlobalElement ( k );
      double_complex z           = computeGl ( globalIndex, *psi );
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
        TEST_FOR_EXCEPTION( true,
    			            std::logic_error,
    			            "Illegal equationType " << eqType );
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
  Teuchos::RCP<ComplexVector > psi = Teuchos::null;

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
        TEST_FOR_EXCEPTION( true,
    			            std::logic_error,
    			            "Illegal equationType " << eqType );
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
  for ( unsigned int k=0; k<psi.getLocalLength(); k++ )
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
double
GinzburgLandau::normalizedScaledL2Norm ( const ComplexVector &psi
                                       ) const
{
  double localSum = 0.0;

  Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

  // sum up on each processor
  for ( unsigned int k=0; k<psi.getLocalLength(); k++ )
    {
	  int kGlobal = psi.getMap()->getGlobalElement ( k );
	  double area = grid_->cellArea( kGlobal );
	  localSum += area * norm ( psiView[k] );
    }

  // reduce and scatter such that energy is available on
  // all cores
  int count = 1; // send *one* integer
  Teuchos::Array<double> sendBuff ( count ), recvBuff ( count );

  // fill send buffer
  sendBuff[0] = localSum;

  int numProcs =  psi.getMap()->getComm()->getSize();
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
  double domainArea = grid_->getGridDomainArea();
  double l2norm = sqrt(recvBuff[0]) / domainArea;

  return l2norm;
}
// =============================================================================
// Count the number of vortices by the total phase change along the boundary
// of the domain.
// TODO Make this work in multicore environments.
int GinzburgLandau::getVorticity ( const ComplexVector &psi
                                 ) const
{
  // this function only works
  int numProcs = psi.getMap()->getComm()->getSize();
  if ( numProcs!=1 )
    return -1;

  int vorticity = 0;

  const double PI = 3.14159265358979323846264338327950288419716939937510;
  // Consider jumps in the argument greater than threshold phase jumps.
  const double threshold = 1.5*PI;

  // Get a view of the whole vector.
  // Remember: This only works with one core.
  Teuchos::ArrayRCP<const double_complex> psiView = psi.get1dView();

  int numBoundaryPoints = grid_->getNumBoundaryPoints();

  int l = 0;
  double angle = std::arg( psiView[ grid_->boundaryIndex2globalIndex(l) ] );
  double angle0 = angle;
  double anglePrevious;
  for ( l=1; l<numBoundaryPoints; l++ ){
	  anglePrevious = angle;
	  angle = std::arg( psiView[grid_->boundaryIndex2globalIndex(l)] );
      if ( angle-anglePrevious<-threshold )
        vorticity++;
      else if ( angle-anglePrevious>threshold )
        vorticity--;
  }

  // close the circle
  anglePrevious = angle;
  if ( angle0-anglePrevious<-threshold )
    vorticity++;
  else if ( angle0-anglePrevious>threshold )
    vorticity--;

  return vorticity;
}
// =============================================================================
void
GinzburgLandau::writeStateToFile( const Teuchos::RCP<const ComplexVector> &psi,
                                  Teuchos::ParameterList &params,
                                  const std::string &filePath
                                ) const
{
  int numUnknowns = psi->getGlobalLength();
  std::vector<int> p(numUnknowns);
  grid_->lexicographic2grid(&p);

  // Create multivector containing the components that we would like to print.
  // Also make sure the entries appear in lexicographic order.
  Teuchos::RCP<const Tpetra::Map<int> > psiMap = psi->getMap();
  Tpetra::MultiVector<double, int> psiSplit(psiMap, 2, true);
  Teuchos::ArrayRCP<const double_complex> psiView = psi->get1dView();
  int globalLength = psiSplit.getGlobalLength();
  for (int k = 0; k < globalLength; k++) {
    psiSplit.replaceLocalValue(k, 0, norm(psiView[p[k]]));
    psiSplit.replaceLocalValue(k, 1, arg(psiView[p[k]]));
  }

  grid_->writeWithGrid( psiSplit, params, filePath );
}
// =============================================================================
void
GinzburgLandau::writeSolutionToFile( const Teuchos::RCP<const ComplexVector> &psi,
                                     const std::string &filePath
                                   ) const
{
  // create a parameter list that contains useful items for a solution file
  Teuchos::ParameterList params;

  // TODO avoid calculating free energy and vorticity twice
  double energy = freeEnergy( *psi );
  int vorticity = getVorticity( *psi );

  params.get("H0", A_->getH0() );
  params.get("edge length", grid_->getEdgeLength());
  params.get("Nx", grid_->getNx());
  params.get("free energy", energy);
  params.get("vorticity", vorticity);

  writeStateToFile( psi, params, filePath );
}
// =============================================================================
void
GinzburgLandau::writeAbstractStateToFile( const Teuchos::RCP<const ComplexVector> &psi,
                                          const std::string &filePath
                                        ) const
{
	Teuchos::ParameterList params;
	params.get("edge length", grid_->getEdgeLength());
	params.get("Nx", grid_->getNx());

	writeStateToFile( psi, params, filePath );
}
// =============================================================================
void
GinzburgLandau::appendStats( std::ofstream & fileStream,
		                     const bool header,
		                     const Teuchos::RCP<const ComplexVector> & psi
		                   ) const
{
    if ( header ) {
    	fileStream << "H0                \t"
    			   << "edge length       \t"
    	 		   << "free energy       \t"
    	           << "||x||_2 scaled    \t"
    	           << "vorticity";
    }
    else {
    	// TODO avoid calculating the free energy twice
    	double energy = freeEnergy( *psi );
    	double l2norm = normalizedScaledL2Norm( *psi );
    	int vorticity = getVorticity( *psi );
    	fileStream << A_->getH0() << " \t"
    			   << grid_->getEdgeLength() << " \t"
    			   << energy << " \t"
    	    	   << l2norm << " \t"
    			   << vorticity;
    }

}
// =============================================================================
// NOT A MEMBER OF GinzburgLandau!
void
readStateFromFile ( const Teuchos::RCP<const Teuchos::Comm<int> > & Comm,
		    const std::string                             & filePath,
                    Teuchos::RCP<ComplexVector>                   & psi,
                    Teuchos::RCP<Grid>                            & grid,
                    Teuchos::ParameterList                        & params
                  )
{

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read the raw data, parameters, grid
  Teuchos::RCP<Tpetra::MultiVector<double, int> > psiLexicographicSplit;

  readWithGrid( Comm, filePath, psiLexicographicSplit, grid, params );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Reorder the input guess and put abs and arg part together again.
  // If there was is an initial guess, make sure to get the ordering correct.
  // TODO Look into having this done by Trilinos. If executed on a multiproc
  //      environment, we don't want p to be fully present on all processors.
  int NumComplexUnknowns = psiLexicographicSplit->getGlobalLength();
  // fill p
  std::vector<int> p(NumComplexUnknowns);
  grid->lexicographic2grid(&p);

  psi = Teuchos::rcp(new ComplexVector(psiLexicographicSplit->getMap()));

  // TODO Replace this with get1dView or get2dView of the full MultiVector
  // TODO Make the following work on multiproc
  Teuchos::ArrayRCP<const double> psiSplitAbsView = psiLexicographicSplit->getVector(0)->get1dView();
  Teuchos::ArrayRCP<const double> psiSplitArgView = psiLexicographicSplit->getVector(1)->get1dView();

  for (int k = 0; k < NumComplexUnknowns; k++) {
      double_complex value = std::polar( sqrt(psiSplitAbsView[k]),
                                         psiSplitArgView[k]
			                             );
       psi->replaceGlobalValue(p[k], value);
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
// =============================================================================
