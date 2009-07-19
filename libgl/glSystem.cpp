#include "glSystem.h"
#include <iostream>
#include <complex>

#include "Epetra_Export.h"

// =============================================================================
// Class constructor
GlSystem::GlSystem( int nx,
                    GinzburgLandau::GinzburgLandau gl,
                    Epetra_Comm& comm ):
  NumGlobalElements(0),
  NumMyElements(0),  // gets set after map creation
  NumComplexUnknowns(0),
  Gl(gl),
  Nx(nx),
  Comm(&comm),
  rhs(0),
  StandardMap(0),
  EverywhereMap(0)
{
  int NumComplexUnknowns = (Nx+1)*(Nx+1);
  int NumGlobalElements = 2*NumComplexUnknowns+1;

  // define the map where the data is nicely distributed over all processors
  StandardMap = new Epetra_Map( NumGlobalElements, 0, *Comm );
  NumMyElements = StandardMap->NumMyElements();

  // define the map where each processor has a full solution vector
  EverywhereMap = new Epetra_Map( NumGlobalElements, NumGlobalElements, 0, *Comm );
}
// =============================================================================


// =============================================================================
// Destructor
GlSystem::~GlSystem()
{
// delete graph
  delete StandardMap;
  delete EverywhereMap;
}
// =============================================================================


// // =============================================================================
// int psiIndex2realIndex ( int psiIndex, complexPart cPart )
// {
// 
//   select (cPart) {
//       case REAL:
//           realIndex = 2*psiIndex-1;
//       case IMAGINARY:
//           realIndex = 2*psiIndex;
//   }
// 
//   return realIndex;
// }
// // =============================================================================


// =============================================================================
int GlSystem::realIndex2psiIndex ( int realIndex )
{
  return realIndex/2;
}
// =============================================================================


// =============================================================================
bool GlSystem::real2psi( Epetra_Vector realvec,
                         std::complex<double>* psi )
{
  for (int k=0; k<NumComplexUnknowns; k++ )
      psi[k] = ( realvec[2*k-1], realvec[2*k] );

  return true;
}
// =============================================================================



// =============================================================================
bool GlSystem::computeF( const Epetra_Vector& x,
                         Epetra_Vector& FVec )
{
  Epetra_Vector xEverywhere ( *EverywhereMap );
  std::complex<double> psi[NumComplexUnknowns];
  std::complex<double> val;

  // scatter x over all processors
  Epetra_Export Exporter( *StandardMap, *EverywhereMap );
  xEverywhere.Export( x, Exporter, Insert );
  real2psi( xEverywhere, psi );

  // loop over the system rows
  for( int i=0; i<NumMyElements; i+2 ) {
      // always extract TWO equations at once (real+imaginary part)
      int myGlobalIndex = StandardMap->GID(i);
      // get the index of the complex valued equation
      int psiIndex = realIndex2psiIndex( myGlobalIndex );
      val = Gl.computeGl( psiIndex, psi );

      // fill in the values
      double passVal = real(val);
      rhs->ReplaceGlobalValues( 1, &passVal, &myGlobalIndex );
      if ( i+1<NumMyElements ) {
          passVal = imag(val);
          int passInd = myGlobalIndex+1;
          rhs->ReplaceGlobalValues( 1, &passVal, &passInd );
      }
  }

  return true;
}
// =============================================================================



// // =============================================================================
// computeJacobian( x )
// {
//   complex<double> val;
// 
//   psi = real2psi(x);
// 
//   // distribute psi (x) to each processor
//   Gl.computeJacobianBlocks( myGlobalRowIndex,
//                             psi,
//                             int* columnIndicesPsi, 
//                             int* columnIndicesPsiConj,
//                             std::complex<double>* valuesPsi,
//                             std::complex<double>* valuesPsiConj )
// 
//   ReplaceValue( myGlobalRowIndex, columnIndicesPsi, valuesPsi );
// 
// }
// // =============================================================================

// // =============================================================================
// createGraph()
// // =============================================================================