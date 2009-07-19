#include "glSystem.h"
#include <iostream>
#include <complex>
#include <vector>

#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"

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
void GlSystem::real2psi( Epetra_Vector realvec,
                         std::complex<double>* psi )
{
  for (int k=0; k<NumComplexUnknowns; k++ )
      psi[k] = ( realvec[2*k-1], realvec[2*k] );
}
// =============================================================================


// =============================================================================
bool GlSystem::computeF( const Epetra_Vector& x,
                         Epetra_Vector& FVec     )
{
  Epetra_Vector xEverywhere ( *EverywhereMap );
  std::complex<double>* psi = new std::complex<double>[NumComplexUnknowns];
  std::complex<double> val;

  // scatter x over all processors
  Epetra_Export Exporter( *StandardMap, *EverywhereMap );
  xEverywhere.Export( x, Exporter, Insert );
  (void) real2psi( xEverywhere, psi );

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



// =============================================================================
bool GlSystem::computeJacobian( const Epetra_Vector& x )
{
  int ierr,
      numEntriesPsi,
      numEntriesPsiConj;
  Epetra_Vector xEverywhere ( *EverywhereMap );
  std::complex<double>* psi = new std::complex<double>[NumComplexUnknowns];

  int *columnIndicesPsi     = NULL,
      *columnIndicesPsiReal = NULL,
      *columnIndicesPsiImag = NULL,
      *columnIndicesPsiConj = NULL,
      *columnIndicesPsiConjReal = NULL,
      *columnIndicesPsiConjImag = NULL;
  std::complex<double> *valuesPsi     = NULL,
                       *valuesPsiConj = NULL;
  double *valuesPsiReal     = NULL,
         *valuesPsiConjReal = NULL,
         *valuesPsiImag     = NULL,
         *valuesPsiConjImag = NULL;

  // scatter x over all processors
  Epetra_Export Exporter( *StandardMap, *EverywhereMap );
  xEverywhere.Export( x, Exporter, Insert );
  real2psi( xEverywhere, psi );

  // Construct the Epetra Matrix
  for( int i=0 ; i<NumMyElements ; i+2 ) {
      int Row = StandardMap->GID(i);
      // get the values and column indices
      Gl.computeJacobianBlocks( Row,
                                psi,
                                numEntriesPsi,
                                columnIndicesPsi,
                                valuesPsi,
                                numEntriesPsiConj,
                                columnIndicesPsiConj,
                                valuesPsiConj );

      // insert the coefficients Re(alpha_i) of Re(psi_i)
      columnIndicesPsiReal = new int [numEntriesPsi];
      for (int k=0; k<numEntriesPsi; k++)
          columnIndicesPsiReal[k] = 2*columnIndicesPsi[k]-1;
      valuesPsiReal = new double(numEntriesPsi);
      for (int k=0; k<numEntriesPsi; k++)
          valuesPsiReal[k] = std::real( valuesPsi[k] );
      ierr = jacobian->ReplaceGlobalValues( Row,
                                            numEntriesPsi,
                                            valuesPsiReal,
                                            columnIndicesPsiReal );

      // insert the coefficients Re(beta_i) of Re(psi_i)
      columnIndicesPsiConjReal = new int[numEntriesPsiConj];
      for (int k=0; k<numEntriesPsiConj; k++)
          columnIndicesPsiConjReal[k] = 2*columnIndicesPsiConj[k]-1;
      valuesPsiConjReal = new double[numEntriesPsiConj];
      for (int k=0; k<numEntriesPsiConj; k++)
          valuesPsiConjReal[k] = std::real( valuesPsiConj[k] );
      ierr = jacobian->SumIntoMyValues( Row,
                                        numEntriesPsiConj,
                                        valuesPsiConjReal,
                                        columnIndicesPsiConjReal );

      // insert the coefficients -Im(alpha_i) of Im(psi_i)
      columnIndicesPsiImag = new int[numEntriesPsi];
      for (int k=0; k<numEntriesPsi; k++)
          columnIndicesPsiImag[k] = 2*columnIndicesPsi[k];
      valuesPsiImag = new double[numEntriesPsi];
      for (int k=0; k<numEntriesPsi; k++)
          valuesPsiImag[k] = -std::imag( valuesPsi[k] );
      ierr = jacobian->ReplaceGlobalValues( Row,
                                            numEntriesPsi,
                                            valuesPsiImag,
                                            columnIndicesPsiImag );

      // insert the coefficients Im(beta_i) of Im(psi_i)
      columnIndicesPsiConjImag = new int[numEntriesPsiConj];
      for (int k=0; k<numEntriesPsiConj; k++)
          columnIndicesPsiConjImag[k] = 2*columnIndicesPsiConj[k];
      valuesPsiConjImag = new double[numEntriesPsiConj];
      for (int k=0; k<numEntriesPsiConj; k++)
          valuesPsiConjImag[k] = std::imag( valuesPsiConj[k] );
      ierr = jacobian->SumIntoMyValues( Row,
                                        numEntriesPsiConj,
                                        valuesPsiConjImag,
                                        columnIndicesPsiConjImag );

//       // -----------------------------------------------------------------
//       if ( i+1<NumMyElements ) {
//       // insert the coefficients Im(alpha_i) of Re(psi_i)
//       int* columnIndicesPsiImag = 2*columnIndicesPsi-1;
//       numEntries = columnIndicesPsiImag.length();
//       double* valuesPsiImag = imag(valuesPsi);
//       ierr = jacobian->ReplaceGlobalValues( Row,
//                                             numEntries,
//                                             &valuesPsiImag,
//                                             &columnIndicesPsiImag );
// 
//       // insert the coefficients Imag(beta_i) of Re(psi_i)
//       int* columnIndicesPsiImag = 2*columnIndicesPsiConj-1;
//       numEntries = columnIndicesPsiConjImag.length();
//       double* valuesPsiConjImag = real(valuesPsiConj);
//       ierr = jacobian->AddToGlobalValues( Row,
//                                           numEntries,
//                                           &valuesPsiConjImag,
//                                           &columnIndicesPsiImag );
// 
//       // insert the coefficients -Im(alpha_i) of Im(psi_i)
//       int* columnIndicesPsiImag = 2*columnIndicesPsi;
//       numEntries = columnIndicesPsiImag.length();
//       double* valuesPsiImag = imag(valuesPsi);
//       ierr = jacobian->ReplaceGlobalValues( Row,
//                                             numEntries,
//                                             &valuesPsiImag,
//                                             &columnIndicesPsiImag );
// 
//       // insert the coefficients Im(beta_i) of Im(psi_i)
//       int* columnIndicesPsiConjImag = 2*columnIndicesPsi;
//       numEntries = columnIndicesPsiConjImag.length();
//       double* valuesPsiConjImag = imag(valuesPsi);
//       ierr = jacobian->AddToGlobalValues( Row,
//                                           numEntries,
//                                           &valuesPsiConjImag,
//                                           &columnIndicesPsiConjImag );
//       }
      // -----------------------------------------------------------------

  }

  return true;
}
// =============================================================================

// // =============================================================================
// createGraph()
// // =============================================================================