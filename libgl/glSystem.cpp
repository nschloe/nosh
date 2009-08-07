#include "glSystem.h"
#include "ioFactory.h"
#include "glException.h"

#include <iostream>
#include <complex>
#include <vector>

#include <Epetra_Export.h>
#include <Epetra_CrsMatrix.h>

#include <EpetraExt_RowMatrixOut.h>

#include <EpetraExt_Utils.h>

// abbreviate the complex type name
typedef std::complex<double> double_complex;

// =============================================================================
// Default constructor
GlSystem::GlSystem( GinzburgLandau::GinzburgLandau &gl,
                    Epetra_Comm                    &comm,
                    std::vector<double_complex>    *psi  ):
  NumGlobalElements(0),
  NumMyElements(0),
  NumComplexUnknowns(0),
  Gl(gl),
  Comm(&comm),
  StandardMap(0),
  EverywhereMap(0),
  rhs(0),
  Graph(0),
  jacobian(0),
  initialSolution(0),
  outputDir(".")
{
  int Nx = Gl.getStaggeredGrid()->getNx();
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

  if ( !psi ) { // null pointer
      initialSolution->PutScalar(0.5); // Default initialization
  } else {
      if ( psi->size() != NumComplexUnknowns ) {
          std::string message = "Size of the initial guess vector ("
                              + EpetraExt::toString( int(psi->size()) )
                              + ") does not coincide with the number of unknowns ("
                              + EpetraExt::toString( NumComplexUnknowns ) + ").";
          throw glException( "GlSystem::GlSystem",
                             message );
      }
      complex2real( *psi, initialSolution );
  }


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
  delete Graph;
  delete StandardMap;
  delete EverywhereMap;
}
// =============================================================================


// =============================================================================
int GlSystem::realIndex2complexIndex ( const int realIndex )
{
  if ( !(realIndex%2) ) // realIndex even
      return realIndex/2;
  else
      return (realIndex-1)/2;
}
// =============================================================================


// =============================================================================
void GlSystem::real2complex( const Epetra_Vector    &realvec,
                             vector<double_complex> &psi      )
{
  for (int k=0; k<NumComplexUnknowns; k++ )
      psi[k] = double_complex( realvec[2*k], realvec[2*k+1] );
}
// =============================================================================



// =============================================================================
void GlSystem::complex2real( const vector<double_complex> &psi,
                             Teuchos::RCP<Epetra_Vector>  realvec )
{
  for (int k=0; k<NumComplexUnknowns; k++ ) {
      (*realvec)[2*k]   = real(psi[k]);
      (*realvec)[2*k+1] = imag(psi[k]);
  }
}
// =============================================================================



// =============================================================================
bool GlSystem::computeF( const Epetra_Vector &x,
                         Epetra_Vector       &FVec,
                         const NOX::Epetra::Interface::Required::FillType fillFlag  )
{
  Epetra_Vector          xEverywhere ( *EverywhereMap );
  vector<double_complex> psi(NumComplexUnknowns);
  double_complex         val;

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
  // compute the values of the Jacobian
  createJacobian( VALUES, x );

  // optimize storage
  jacobian->FillComplete();

  // Sync up processors to be safe
  Comm->Barrier();

  return true;
}
// =============================================================================


// ==========================================================================
bool GlSystem::computePreconditioner( const Epetra_Vector    &x,
                                      Epetra_Operator        &Prec,
                                      Teuchos::ParameterList *precParams )
{
  std::cerr << "ERROR: GlSystem::preconditionVector() - "
            << "    Use Explicit Jacobian only for this test problem!" << std::endl;
  throw "GlSystem Error";
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
  std::vector<double_complex> psi(NumComplexUnknowns);

  vector<int>            colIndA, colIndB;
  vector<double_complex> valuesA, valuesB;

  int    *colInd      = NULL,
         *colIndAReal = NULL, *colIndAImag = NULL,
         *colIndBReal = NULL, *colIndBImag = NULL;
  double *values      = NULL,
         *valuesAReal = NULL, *valuesAImag = NULL,
         *valuesBReal = NULL, *valuesBImag = NULL;

  int ierr, k, complexRow, numEntries;

  if (jc==VALUES) {
      Epetra_Vector xEverywhere ( *EverywhereMap );
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
          colInd = new int[numEntries];
          for (int k=0; k<numEntries; k++ )
              colInd[k] = k;
          if (jc==VALUES) {// fill on columns and values
              values = new double[numEntries];
              for (int k=0; k<NumComplexUnknowns; k++ ) {
                  values[2*k]   = -imag(psi[k]);
                  values[2*k+1] =  real(psi[k]);
              }
              // fill it in!
              ierr = jacobian->SumIntoGlobalValues( Row,
                                                    numEntries,
                                                    values,
                                                    colInd );
              delete [] values;
              values = NULL;
          } else {// only fill the sparsity graph
              Graph->InsertGlobalIndices( Row,
                                          numEntries,
                                          colInd );
          }
          delete [] colInd;
          colInd = NULL;


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
                                 colIndA, valuesA,
                                 colIndB, valuesB );
          } else {
              // only fill the sparsity graph
              Gl.getJacobianRowSparsity( complexRow,
                                         colIndA,
                                         colIndB );
          }

          if ( !(Row%2) ) { // myGlobalIndex is even <=> Real part of the equation system
              // ---------------------------------------------------------------
              // insert the coefficients Re(A) of Re(psi)
              numEntries = colIndA.size();
              colIndAReal = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  colIndAReal[k] = 2*colIndA[k];

              if (jc==VALUES) {
                  valuesAReal = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesAReal[k] = real( valuesA[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesAReal,
                                                        colIndAReal );
                  delete [] valuesAReal; valuesAReal = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              colIndAReal );
              }
              delete [] colIndAReal; colIndAReal = NULL;

              // insert the coefficients Re(B) of Re(psi)
              numEntries = colIndB.size();
              colIndBReal = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  colIndBReal[k] = 2*colIndB[k];
              if (jc==VALUES) {
                  valuesBReal = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesBReal[k] = real( valuesB[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesBReal,
                                                        colIndBReal );
                  delete [] valuesBReal;
                  valuesBReal = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              colIndBReal );
              }
              delete [] colIndBReal;
              colIndBReal = NULL;

              // insert the coefficients -Im(A) of Im(psi)
              numEntries = colIndA.size();
              colIndAImag = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  colIndAImag[k] = 2*colIndA[k]+1;
              if (jc==VALUES) {
                  valuesAImag = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesAImag[k] = -imag( valuesA[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesAImag,
                                                        colIndAImag );
                  delete [] valuesAImag;
                  valuesAImag = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              colIndAImag );
              }
              delete [] colIndAImag;
              colIndAImag = NULL;


              // insert the coefficients Im(B) of Im(psi)
              numEntries = colIndB.size();
              colIndBImag = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  colIndBImag[k] = 2*colIndB[k]+1;
              if (jc==VALUES) {
                  valuesBImag = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesBImag[k] = imag( valuesB[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesBImag,
                                                        colIndBImag );
                  delete [] valuesBImag;
                  valuesBImag = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              colIndBImag );
              }
              delete [] colIndBImag;
              colIndBImag = NULL;


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
              // insert the coefficients Im(A) of Re(psi)
              numEntries = colIndA.size();
              colIndAReal = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  colIndAReal[k] = 2*colIndA[k];
              if (jc==VALUES) {
                  valuesAImag = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesAImag[k] = imag( valuesA[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesAImag,
                                                        colIndAReal );
                  delete [] valuesAImag;
                  valuesAImag = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              colIndAReal );
              }
              delete [] colIndAReal;
              colIndAReal = NULL;

              // insert the coefficients Im(B) of Re(psi)
              numEntries = colIndB.size();
              colIndBReal = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  colIndBReal[k] = 2*colIndB[k];
              if (jc==VALUES) {
                  valuesBImag = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesBImag[k] = imag( valuesB[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesBImag,
                                                        colIndBReal );
                  delete [] valuesBImag;
                  valuesBImag = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              colIndBReal );
              }
              delete [] colIndBReal;
              colIndBReal = NULL;


              // insert the coefficients Re(A) of Im(psi)
              numEntries = colIndA.size();
              colIndAImag = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  colIndAImag[k] = 2*colIndA[k]+1;
              if (jc==VALUES) {
                  valuesAReal = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesAReal[k] = real( valuesA[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesAReal,
                                                        colIndAImag );
                  delete [] valuesAReal;
                  valuesAReal = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              colIndAImag );
              }
              delete [] colIndAImag;
              colIndAImag = NULL;


              // insert the coefficients -Re(B) of Im(psi)
              numEntries = colIndB.size();
              colIndBImag = new int[numEntries];
              for (k=0; k<numEntries; k++)
                  colIndBImag[k] = 2*colIndB[k]+1;
              if (jc==VALUES) {
                  valuesBReal = new double[numEntries];
                  for (k=0; k<numEntries; k++)
                      valuesBReal[k] = -real( valuesB[k] );
                  ierr = jacobian->SumIntoGlobalValues( Row,
                                                        numEntries,
                                                        valuesBReal,
                                                        colIndBImag );
                  delete [] valuesBReal;
                  valuesBReal = NULL;
              } else {
                  Graph->InsertGlobalIndices( Row,
                                              numEntries,
                                              colIndBImag );
              }
              delete [] colIndBImag; colIndBImag = NULL;

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
          std::string message = "FillComplete returned error code "
                              + EpetraExt::toString( i ) + ".";
          throw glException( "GlSystem::createGraph",
                             message );
  }
  // ---------------------------------------------------------------------------

  // Sync up processors for safety's sake
  Comm->Barrier();

  return true;
}
// =============================================================================


// =============================================================================
// function used by LOCA
void GlSystem::setParameters(const LOCA::ParameterVector &p)
{
  double h0 = p.getValue("H0");

  // set H0 in the underlying problem class
  Gl.getStaggeredGrid()->setH0( h0 );
}
// =============================================================================




// ===========================================================================
void GlSystem::dataForPrintSolution( const int conStep_,
                                     const int timeStep_,
                                     const int totalTimeSteps_ )
{
//     myConStep = conStep_;
    std::cout << "ggggggggggg" << std::endl;
}
// ===========================================================================




// =============================================================================
// function used by LOCA
void GlSystem::printSolution( const Epetra_Vector &x,
                              double              conParam )
{
  // convert the real valued vector to psi
  vector<double_complex> psi(NumComplexUnknowns);
  real2complex( x, psi );

  // create temporary parameter list
  // TODO: get rid of this
  Teuchos::ParameterList tmpList;
  tmpList.get( "ContinuatonParameter", conParam );

  std::string fileName = "cont.vtk";

  // TODO: print other parameters as well
  IoVirtual* fileIo = IoFactory::createFileIo( outputDir+"/"+fileName );
  fileIo->write( psi,
                 tmpList,
                 *(Gl.getStaggeredGrid()) );
}
// =============================================================================



// =============================================================================
// function used by LOCA
void GlSystem::setOutputDir( const string &directory )
{
   outputDir = directory;
}
// =============================================================================




// =============================================================================
void GlSystem::solutionToFile( const Epetra_Vector          &x,
                               const Teuchos::ParameterList &problemParams,
                               const std::string            &fileName       )
{
  // convert the real valued vector to psi
  vector<double_complex> psi(NumComplexUnknowns);
  real2complex( x, psi );

  IoVirtual* fileIo = IoFactory::createFileIo( fileName );
  fileIo->write( psi,
                 problemParams,
                 *(Gl.getStaggeredGrid()) );

}
// =============================================================================