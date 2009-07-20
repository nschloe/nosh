#include "glSystem.h"
#include <iostream>
#include <complex>
#include <vector>

#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"

// =============================================================================
// Class constructor
GlSystem::GlSystem( int nx,
                    double h0,
                    double edgelength,
                    Epetra_Comm& comm ):
  NumGlobalElements(0),
  NumMyElements(0),  // gets set after map creation
  NumComplexUnknowns(0),
  Nx(nx),
  H0(h0),
  Gl(GinzburgLandau::GinzburgLandau( nx, edgelength, h0 )),
  Edgelength(edgelength),
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

  // initialize solution
  initialSolution = Teuchos::rcp(new Epetra_Vector(*StandardMap));
  initializeSoln();

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
                         Epetra_Vector& FVec,
                         const NOX::Epetra::Interface::Required::FillType fillFlag  )
{
  Epetra_Vector xEverywhere ( *EverywhereMap );
  std::complex<double>* psi = new std::complex<double>[NumComplexUnknowns];
  std::complex<double> val;

  // scatter x over all processors
  Epetra_Export Exporter( *StandardMap, *EverywhereMap );
  xEverywhere.Export( x, Exporter, Insert );
  (void) real2psi( xEverywhere, psi );

  // loop over the system rows
  double passVal;
  for( int i=0; i<NumMyElements; i++ ) {
      int myGlobalIndex = StandardMap->GID(i);
      if (myGlobalIndex==2*NumComplexUnknowns+1) { // phase condition
          passVal = 0.0;
      } else { // GL equations
          // get the index of the complex valued equation
          int psiIndex = realIndex2psiIndex( myGlobalIndex );
          // get the complex value
          // TODO: The same value is actually fetched twice in this loop
          //       possibly consecutively: Once for the real, once for
          //       the imaginary part of it.
          val = Gl.computeGl( psiIndex, psi );

          if ( myGlobalIndex%2 ) // myGlobalIndex is even
              passVal = real(val);
          else // myGlobalIndex is odd
              passVal = imag(val);
      }
      rhs->ReplaceGlobalValues( 1, &passVal, &myGlobalIndex );
  }

  return true;
}
// =============================================================================



// =============================================================================
bool GlSystem::computeJacobian( const Epetra_Vector& x,
                                Epetra_Operator &Jac    )
{
  int ierr,
      numEntriesPsi,
      numEntriesPsiConj;
  Epetra_Vector xEverywhere ( *EverywhereMap );
  std::complex<double>* psi = new std::complex<double>[NumComplexUnknowns];

  int *columnIndices        = NULL,
      *columnIndicesPsi     = NULL,
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
  for( int i=0 ; i<NumMyElements ; i++ ) {
      int Row = StandardMap->GID(i);

      if (Row==2*NumComplexUnknowns+1) {// phase condition
          // fill in phase condition stuff
          int numEntries = 2*NumComplexUnknowns;
          columnIndices = new int[numEntries];
          for (int k; k<numEntries; k++ )
              columnIndices[k] = k;
          double* values = new double[numEntries];
          for (int k; k<NumComplexUnknowns; k++ ) {
              values[2*k]   = -imag(psi[k]);
              values[2*k+1] =  real(psi[k]);
          }
          // fill it in!
          ierr = jacobian->ReplaceGlobalValues( Row,
                                                numEntriesPsi,
                                                values,
                                                columnIndices );
      } else { // GL equations
          // get the values and column indices
          // TODO: The same value is actually fetched twice in this loop
          //       possibly consecutively: Once for the real, once for
          //       the imaginary part of it.
          Gl.computeJacobianBlocks( Row,
                                    psi,
                                    numEntriesPsi,
                                    columnIndicesPsi,
                                    valuesPsi,
                                    numEntriesPsiConj,
                                    columnIndicesPsiConj,
                                    valuesPsiConj );

          if ( Row%2 ) { // myGlobalIndex is even
              // ---------------------------------------------------------
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

              // right bordering
              int k = realIndex2psiIndex( Row );
              double value  = imag(psi[k]);
              int    column = 2*NumComplexUnknowns+1;
              ierr = jacobian->ReplaceMyValues( Row,
                                                1,
                                                &value,
                                                &column );

              // ---------------------------------------------------------
          } else { // myGlobalIndex is odd
              // ---------------------------------------------------------
              // insert the coefficients Im(alpha_i) of Re(psi_i)
              columnIndicesPsiReal = new int [numEntriesPsi];
              for (int k=0; k<numEntriesPsi; k++)
                  columnIndicesPsiReal[k] = 2*columnIndicesPsi[k]-1;
              valuesPsiImag = new double(numEntriesPsi);
              for (int k=0; k<numEntriesPsi; k++)
                  valuesPsiImag[k] = std::imag( valuesPsi[k] );
              ierr = jacobian->ReplaceGlobalValues( Row,
                                                    numEntriesPsi,
                                                    valuesPsiImag,
                                                    columnIndicesPsiReal );

              // insert the coefficients Im(beta_i) of Re(psi_i)
              columnIndicesPsiConjReal = new int[numEntriesPsiConj];
              for (int k=0; k<numEntriesPsiConj; k++)
                  columnIndicesPsiConjReal[k] = 2*columnIndicesPsiConj[k]-1;
              valuesPsiConjImag = new double[numEntriesPsiConj];
              for (int k=0; k<numEntriesPsiConj; k++)
                  valuesPsiConjImag[k] = std::imag( valuesPsiConj[k] );
              ierr = jacobian->SumIntoMyValues( Row,
                                                numEntriesPsiConj,
                                                valuesPsiConjImag,
                                                columnIndicesPsiConjReal );

              // insert the coefficients Re(alpha_i) of Im(psi_i)
              columnIndicesPsiImag = new int[numEntriesPsi];
              for (int k=0; k<numEntriesPsi; k++)
                  columnIndicesPsiImag[k] = 2*columnIndicesPsi[k];
              valuesPsiReal = new double[numEntriesPsi];
              for (int k=0; k<numEntriesPsi; k++)
                  valuesPsiReal[k] = std::imag( valuesPsi[k] );
              ierr = jacobian->ReplaceGlobalValues( Row,
                                                    numEntriesPsi,
                                                    valuesPsiReal,
                                                    columnIndicesPsiImag );

              // insert the coefficients -Re(beta_i) of Im(psi_i)
              columnIndicesPsiConjImag = new int[numEntriesPsiConj];
              for (int k=0; k<numEntriesPsiConj; k++)
                  columnIndicesPsiConjImag[k] = 2*columnIndicesPsiConj[k];
              valuesPsiConjReal = new double[numEntriesPsiConj];
              for (int k=0; k<numEntriesPsiConj; k++)
                  valuesPsiConjReal[k] = -std::imag( valuesPsiConj[k] );
              ierr = jacobian->SumIntoMyValues( Row,
                                                numEntriesPsiConj,
                                                valuesPsiConjReal,
                                                columnIndicesPsiConjImag );

              // right bordering
              int k = realIndex2psiIndex( Row );
              double value  = -real(psi[k]);
              int    column = 2*NumComplexUnknowns+1;
              ierr = jacobian->ReplaceMyValues( Row,
                                                1,
                                                &value,
                                                &column );
              // ---------------------------------------------------------
          }
      }
  }

  return true;
}
// =============================================================================


// ==========================================================================
bool GlSystem::computePreconditioner( const Epetra_Vector& x,
                                      Epetra_Operator& Prec,
                                      Teuchos::ParameterList* precParams )
{
  cout << "ERROR: GlSystem::preconditionVector() - "
       << "Use Explicit Jacobian only for this test problem!" << endl;
  throw "GlSystem Error";
}
// ==========================================================================


// ==========================================================================
// Set initialSolution to desired initial condition
bool GlSystem::initializeSoln()
{
  initialSolution->PutScalar(0.0); // Default initialization

//   int n=(NumGlobalElements-1)/2;
//   for (int k=0;k<n;k++){
//     (*initialSolution)[k] = k+1;
//   }
//   for (int k=n;k<2*n;k++){
//     (*initialSolution)[k] = k+1;
//   }
//
//   (*initialSolution)[NumGlobalElements-1] = 2.718291828;

  return true;
}
// ==========================================================================


// ==========================================================================
Teuchos::RCP<Epetra_Vector> GlSystem::getSolution()
{
  return initialSolution;
}
// ==========================================================================


// ==========================================================================
Teuchos::RCP<Epetra_CrsMatrix> GlSystem::getJacobian()
{
  return jacobian;
}
// ==========================================================================



// // =============================================================================
// createGraph()
// // =============================================================================