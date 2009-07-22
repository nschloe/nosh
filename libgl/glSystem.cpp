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
  Gl(GinzburgLandau::GinzburgLandau( nx, edgelength, h0 )),
  Nx(nx),
  H0(h0),
  Edgelength(edgelength),
  Comm(&comm),
  StandardMap(0),
  EverywhereMap(0),
  rhs(0),
  Graph(0),
  jacobian(0),
  initialSolution(0)
{
  NumComplexUnknowns = (Nx+1)*(Nx+1);
  NumGlobalElements  = 2*NumComplexUnknowns+1;

  // define the map where the data is nicely distributed over all processors
  StandardMap = new Epetra_Map( NumGlobalElements, 0, *Comm );
  NumMyElements = StandardMap->NumMyElements();

  // define the map where each processor has a full solution vector
  EverywhereMap = new Epetra_Map( NumGlobalElements,
                                  NumGlobalElements,
                                  0,
                                  *Comm              );

  // initialize solution
  initialSolution = Teuchos::rcp(new Epetra_Vector(*StandardMap));
  initializeSoln();

  // create the sparsity structure (graph) of the Jacobian
  // use x as DUMMY argument
  Epetra_Vector dummy(*StandardMap);
  createJacobian( ONLY_GRAPH, dummy );

  // Allocate the sparsity pattern of the jacobian matrix
  jacobian = Teuchos::rcp(new Epetra_CrsMatrix (Copy, *Graph));

  // Clean-up
  jacobian->FillComplete();
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
int GlSystem::realIndex2complexIndex ( const int realIndex )
{
  return realIndex/2;
}
// =============================================================================


// =============================================================================
void GlSystem::real2complex( const Epetra_Vector            realvec,
                             vector<std::complex<double> >& psi      )
{
  for (int k=0; k<NumComplexUnknowns; k++ )
      psi[k] = std::complex<double>( realvec[2*k], realvec[2*k+1] );
}
// =============================================================================


// =============================================================================
bool GlSystem::computeF( const Epetra_Vector &x,
                         Epetra_Vector       &FVec,
                         const NOX::Epetra::Interface::Required::FillType fillFlag  )
{
  Epetra_Vector xEverywhere ( *EverywhereMap );
  vector<std::complex<double> > psi(NumComplexUnknowns);
  std::complex<double> val;

  // scatter x over all processors
  Epetra_Export Exporter( *StandardMap, *EverywhereMap );
  xEverywhere.Export( x, Exporter, Insert );
  (void) real2complex( xEverywhere, psi );

  // loop over the system rows
  double passVal;
  for(int i=0; i<NumMyElements; i++ ) {
      int myGlobalIndex = StandardMap->GID(i);
      if (myGlobalIndex==2*NumComplexUnknowns) { // phase condition
          passVal = 0.0;
      } else { // GL equations
          // get the index of the complex valued equation
          int psiIndex = realIndex2complexIndex( myGlobalIndex );
          // get the complex value
          // TODO: The same value is actually fetched twice in this loop
          //       possibly consecutively: Once for the real, once for
          //       the imaginary part of it.
          val = Gl.computeGl( psiIndex, psi );

          if ( !(myGlobalIndex%2) ) // myGlobalIndex is even
              passVal = real(val);
          else // myGlobalIndex is odd
              passVal = imag(val);
      }
      FVec.ReplaceGlobalValues( 1, &passVal, &myGlobalIndex );
  }

  return true;
}
// =============================================================================



// =============================================================================
bool GlSystem::computeJacobian( const Epetra_Vector &x,
                                Epetra_Operator     &Jac    )
{

// cout << "x:" << endl;
// cout << x << endl;

  // compute the values of the Jacobian
  createJacobian( VALUES, x );

  // optimize storage
  jacobian->FillComplete();

  // Sync up processors to be safe
  Comm->Barrier();

// cout << "jacobian:" << endl;
// cout << *jacobian << endl;

  return true;
}
// =============================================================================


// ==========================================================================
bool GlSystem::computePreconditioner( const Epetra_Vector    &x,
                                      Epetra_Operator        &Prec,
                                      Teuchos::ParameterList *precParams )
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
  initialSolution->PutScalar(0.5); // Default initialization

//   for (int k=0;k<NumComplexUnknowns;k++){
//     (*initialSolution)[2*k]   = k*0.1;
//     (*initialSolution)[2*k+1] = 0.0;
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


// =============================================================================
bool GlSystem::createJacobian( const jacCreator    jc,
                               const Epetra_Vector &x  )
{
  int ierr;
  Epetra_Vector xEverywhere ( *EverywhereMap );
  std::vector<std::complex<double> > psi(NumComplexUnknowns);

  vector<int> columnIndicesPsi,
              columnIndicesPsiConj;

  vector<std::complex<double> > valuesPsi,
                                valuesPsiConj;

  int *columnIndices        = NULL,
      *columnIndicesPsiReal = NULL,
      *columnIndicesPsiImag = NULL,
      *columnIndicesPsiConjReal = NULL,
      *columnIndicesPsiConjImag = NULL;
  double *valuesPsiReal     = NULL,
         *valuesPsiConjReal = NULL,
         *valuesPsiImag     = NULL,
         *valuesPsiConjImag = NULL;

  int k;
  int complexRow, numEntries;

  if (jc==VALUES) {
      // scatter x over all processors
      Epetra_Export Exporter( *StandardMap, *EverywhereMap );
      xEverywhere.Export( x, Exporter, Insert );
      real2complex( xEverywhere, psi );

      // set the matrix to 0
      jacobian->PutScalar(0.0);
  } else {
      if (Graph != 0) {
        delete Graph;
        Graph = 0;
      }

      // allocate the graph
      Graph = new Epetra_CrsGraph( Copy, *StandardMap, 1 );
  }


  // Construct the Epetra Matrix
  for( int i=0 ; i<NumMyElements ; i++ ) {
      int Row = StandardMap->GID(i);

      if (Row==2*NumComplexUnknowns) {// phase condition
          // fill in phase condition stuff
          numEntries = 2*NumComplexUnknowns;
          columnIndices = new int[numEntries];
          for (int k=0; k<numEntries; k++ )
              columnIndices[k] = k;
          if (jc==VALUES) {// fill on columns and values
              double* values = new double[numEntries];
              for (int k=0; k<NumComplexUnknowns; k++ ) {
                  values[2*k]   = -imag(psi[k]);
                  values[2*k+1] =  real(psi[k]);
              }
              // fill it in!
              ierr = jacobian->SumIntoGlobalValues( Row,
                                                    numEntries,
                                                    values,
                                                    columnIndices );
              delete [] values; values = NULL;
          } else {// only fill the sparsity graph
              Graph->InsertGlobalIndices( Row,
                                          numEntries,
                                          columnIndices );
          }
          delete [] columnIndices; columnIndices = NULL;


      } else { // GL equations
          // get the values and column indices
          // TODO: The same value is actually fetched twice in this loop
          //       possibly consecutively: Once for the real, once for
          //       the imaginary part of it.
          if ( !(Row%2) ) // Row even
              complexRow = Row/2;
          else
              complexRow = (Row-1)/2;

          if (jc==VALUES) {
              // fill on columns and values
              Gl.getJacobianRow( complexRow,
                                 psi,
                                 columnIndicesPsi,
                                 valuesPsi,
                                 columnIndicesPsiConj,
                                 valuesPsiConj );
          } else {
              // only fill the sparsity graph
              Gl.getJacobianRowSparsity( complexRow,
                                         columnIndicesPsi,
                                         columnIndicesPsiConj );
          }

          if ( !(Row%2) ) { // myGlobalIndex is even
              // ---------------------------------------------------------
              // insert the coefficients Re(alpha_i) of Re(psi_i)
              numEntries = columnIndicesPsi.size();
              columnIndicesPsiReal = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  columnIndicesPsiReal[k] = 2*columnIndicesPsi[k];

              if (jc==VALUES) {
                  valuesPsiReal = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesPsiReal[k] = std::real( valuesPsi[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesPsiReal,
                                                        columnIndicesPsiReal );
                  delete [] valuesPsiReal; valuesPsiReal = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              columnIndicesPsiReal );
              }
              delete [] columnIndicesPsiReal; columnIndicesPsiReal = NULL;

              // insert the coefficients Re(beta_i) of Re(psi_i)
              numEntries = columnIndicesPsiConj.size();
              columnIndicesPsiConjReal = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  columnIndicesPsiConjReal[k] = 2*columnIndicesPsiConj[k];
              if (jc==VALUES) {
                  valuesPsiConjReal = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesPsiConjReal[k] = real( valuesPsiConj[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesPsiConjReal,
                                                        columnIndicesPsiConjReal );
                  delete [] valuesPsiConjReal;
                  valuesPsiConjReal = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              columnIndicesPsiConjReal );
              }
              delete [] columnIndicesPsiConjReal;
              columnIndicesPsiConjReal = NULL;

              // insert the coefficients -Im(alpha_i) of Im(psi_i)
              numEntries = columnIndicesPsi.size();
              columnIndicesPsiImag = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  columnIndicesPsiImag[k] = 2*columnIndicesPsi[k]+1;
              if (jc==VALUES) {
                  valuesPsiImag = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesPsiImag[k] = -std::imag( valuesPsi[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesPsiImag,
                                                        columnIndicesPsiImag );
                  delete [] valuesPsiImag;
                  valuesPsiImag = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              columnIndicesPsiImag );
              }
              delete [] columnIndicesPsiImag;
              columnIndicesPsiImag = NULL;


              // insert the coefficients Im(beta_i) of Im(psi_i)
              numEntries = columnIndicesPsiConj.size();
              columnIndicesPsiConjImag = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  columnIndicesPsiConjImag[k] = 2*columnIndicesPsiConj[k]+1;
              if (jc==VALUES) {
                  valuesPsiConjImag = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesPsiConjImag[k] = std::imag( valuesPsiConj[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesPsiConjImag,
                                                        columnIndicesPsiConjImag );
                  delete [] valuesPsiConjImag;
                  valuesPsiConjImag = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              columnIndicesPsiConjImag );
              }
              delete [] columnIndicesPsiConjImag;
              columnIndicesPsiConjImag = NULL;


              // right bordering
              int k = realIndex2complexIndex( Row );
              int column = 2*NumComplexUnknowns;
              if (jc==VALUES) {
                  double value  = imag(psi[k]);
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        1,
                                                        &value,
                                                        &column );
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              1,
                                              &column );
              }

              // ---------------------------------------------------------
          } else { // Row is odd <=> Imaginary part of the equation
              // ---------------------------------------------------------
              // insert the coefficients Im(alpha_i) of Re(psi_i)
              numEntries = columnIndicesPsi.size();
              columnIndicesPsiReal = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  columnIndicesPsiReal[k] = 2*columnIndicesPsi[k];
              if (jc==VALUES) {
                  valuesPsiImag = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesPsiImag[k] = imag( valuesPsi[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesPsiImag,
                                                        columnIndicesPsiReal );
                  delete [] valuesPsiImag;
                  valuesPsiImag = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              columnIndicesPsiReal );
              }
              delete [] columnIndicesPsiReal;
              columnIndicesPsiReal = NULL;


              // insert the coefficients Im(beta_i) of Re(psi_i)
              numEntries = columnIndicesPsiConj.size();
              columnIndicesPsiConjReal = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  columnIndicesPsiConjReal[k] = 2*columnIndicesPsiConj[k];
              if (jc==VALUES) {
                  valuesPsiConjImag = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesPsiConjImag[k] = std::imag( valuesPsiConj[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesPsiConjImag,
                                                        columnIndicesPsiConjReal );
                  delete [] valuesPsiConjImag;
                  valuesPsiConjImag = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              columnIndicesPsiConjReal );
              }
              delete [] columnIndicesPsiConjReal;
              columnIndicesPsiConjReal = NULL;


              // insert the coefficients Re(alpha_i) of Im(psi_i)
              numEntries = columnIndicesPsi.size();
              columnIndicesPsiImag = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  columnIndicesPsiImag[k] = 2*columnIndicesPsi[k]+1;
              if (jc==VALUES) {
                  valuesPsiReal = new double[columnIndicesPsi.size()];
                  for (k=0; k<numEntries; k++)
                      valuesPsiReal[k] = std::real( valuesPsi[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesPsiReal,
                                                        columnIndicesPsiImag );
                  delete [] valuesPsiReal;
                  valuesPsiReal = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              columnIndicesPsiImag );
              }
              delete [] columnIndicesPsiImag;
              columnIndicesPsiImag = NULL;


              // insert the coefficients -Re(beta_i) of Im(psi_i)
              numEntries = columnIndicesPsiConj.size();
              columnIndicesPsiConjImag = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  columnIndicesPsiConjImag[k] = 2*columnIndicesPsiConj[k]+1;
              if (jc==VALUES) {
                  valuesPsiConjReal = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesPsiConjReal[k] = -std::real( valuesPsiConj[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesPsiConjReal,
                                                        columnIndicesPsiConjImag );
                  delete [] valuesPsiConjReal;
                  valuesPsiConjReal = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              columnIndicesPsiConjImag );
              }
              delete [] columnIndicesPsiConjImag; columnIndicesPsiConjImag = NULL;

              // right bordering
              int column = 2*NumComplexUnknowns;
              if (jc==VALUES) {
                  int k = realIndex2complexIndex( Row );
                  double value  = -real(psi[k]);
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        1,
                                                        &value,
                                                        &column );
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              1,
                                              &column );
              }
              // ---------------------------------------------------------
          }
      }
  }

  // ---------------------------------------------------------------------------
  // finish up the graph construction
  try{
      Graph->FillComplete();
  }
  catch (int i){
     std::cerr << "createGraph:" << std::endl
               << "    FillComplete returned error code " << i << ". Abort."
               << std::endl;
     exit(EXIT_FAILURE);
  }
  // ---------------------------------------------------------------------------

  // Sync up processors to be safe
  Comm->Barrier();

// cout << "jacobian:" << endl;
// cout << *jacobian << endl;

  return true;
}
// =============================================================================