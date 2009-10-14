#include "glSystem.h"
#include "ioFactory.h"
#include "glException.h"

#include <iostream>
#include <fstream>

#include <complex>
#include <vector>

#include <Epetra_Export.h>
#include <Epetra_CrsMatrix.h>
#include <NOX_Utils.H>

#include <EpetraExt_RowMatrixOut.h>

#include <EpetraExt_Utils.h>

#include <Epetra_Map.h>

#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Thyra_EpetraThyraWrappers.hpp>

#include <Teuchos_DefaultComm.hpp>

// abbreviate the complex type name
typedef std::complex<double> double_complex;

// =============================================================================
// Default constructor
GlSystem::GlSystem ( GinzburgLandau::GinzburgLandau                               &gl,
		     const Teuchos::RCP<Epetra_Comm>                              eComm,
                     const bool                                                   &rv,
                     const Teuchos::RCP<Tpetra::MultiVector<double_complex,int> > psi  ) :
    NumRealUnknowns ( 0 ),
    NumMyElements ( 0 ),
    NumComplexUnknowns ( 0 ),
    Gl ( gl ),
    EComm ( eComm ),
    TComm ( 0 ),
    RealMap ( 0 ),
    ComplexMap ( 0 ),
    rhs ( 0 ),
    Graph ( 0 ),
    jacobian ( 0 ),
    initialSolution ( 0 ),
    reverse ( rv ),
    outputDir ( "." )
{
  int Nx = Gl.getStaggeredGrid()->getNx();
  NumComplexUnknowns = ( Nx+1 ) * ( Nx+1 );
  NumRealUnknowns    = 2*NumComplexUnknowns+1;

  // @TODO
  // There is (until now?) no way to convert a Teuchos::Comm (of psi) 
  // to an Epetra_Comm (of the real valued representation of psi), so the
  // Epetra_Comm has to be generated explicitly, and two communicators are kept
  // side by side all the time. One must make sure that the two are actually
  // equivalent, which can be checked by Thyra's conversion method create_Comm.
  // @TODO
  // Is is actually necessary to have equivalent communicators on the
  // real-valued and the complex-valued side?
  // How to compare two communicators anyway?
  
  // create fitting Tpetra::Comm
  TComm = Thyra::create_Comm( EComm );
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // define maps
  if ( psi.is_null() )   // null pointer
    {
      // define uniform distribution
      ComplexMap = Teuchos::rcp(new Tpetra::Map<int>( NumComplexUnknowns, 0, TComm ) );
    }
  else
    {
      // psi->getMap() returns a CONST map
      ComplexMap = Teuchos::RCP<const Tpetra::Map<int> >( psi->getMap() );
    }

  // get the map for the real values
  makeRealMap( ComplexMap, RealMap );
  
  // set the number of local elements
  NumMyElements = RealMap->NumMyElements();
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Create a map for the real-valued vector to be spread all over all
  // processors.
  // @TODO Remove (the need for) this.
  // define the map where each processor has a full solution vector
  EverywhereMap = Teuchos::rcp( new Epetra_Map ( NumRealUnknowns,
                                                 NumRealUnknowns,
                                                 0,
                                                 *EComm ) );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // initialize solution
  initialSolution = Teuchos::rcp ( new Epetra_Vector ( *RealMap ) );
  if ( psi.is_null() )   // null pointer
    {
      // define map for psi
      initialSolution->PutScalar ( 0.5 ); // Default initialization
    }
  else
    {
      if ( psi->getGlobalLength() != (unsigned int)NumComplexUnknowns )
        {
          std::string message = "Size of the initial guess vector ("
                                + EpetraExt::toString ( int ( psi->getGlobalLength() ) )
                                + ") does not coincide with the number of unknowns ("
                                + EpetraExt::toString ( NumComplexUnknowns ) + ").";
          throw glException ( "GlSystem::GlSystem",
                              message );
        }
      psi2real ( *psi, *initialSolution );
    }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  


  // create the sparsity structure (graph) of the Jacobian
  // use x as DUMMY argument
  Epetra_Vector dummy ( *RealMap );
  createJacobian ( ONLY_GRAPH, dummy );

  // Allocate the sparsity pattern of the jacobian matrix
  jacobian = Teuchos::rcp ( new Epetra_CrsMatrix ( Copy, *Graph ) );

  // Clean-up
  jacobian->FillComplete();
}
// =============================================================================
// Destructor
GlSystem::~GlSystem()
{
}
// =============================================================================
int GlSystem::realIndex2complexIndex ( const int realIndex )
{
  if ( ! ( realIndex%2 ) ) // realIndex even
    return realIndex/2;
  else
    return ( realIndex-1 ) /2;
}
// =============================================================================
void GlSystem::real2complex ( const Epetra_Vector    &realvec,
                              vector<double_complex> &psi )
{
  for ( int k=0; k<NumComplexUnknowns; k++ )
    psi[k] = double_complex ( realvec[2*k], realvec[2*k+1] );
}
// =============================================================================
void GlSystem::complex2real ( const vector<double_complex> &psi,
                              Teuchos::RCP<Epetra_Vector>  realvec )
{
  for ( int k=0; k<NumComplexUnknowns; k++ )
    {
      ( *realvec ) [2*k]   = real ( psi[k] );
      ( *realvec ) [2*k+1] = imag ( psi[k] );
    }
}
// =============================================================================
// converts a real-valued vector to a complex-valued psi vector
void GlSystem::real2psi ( const Epetra_Vector                     &realvec,
                          Tpetra::MultiVector<double_complex,int> &psi     )
{
  for ( unsigned int k=0; k<psi.getGlobalLength(); k++ ) {
    double_complex z = double_complex ( realvec[2*k], realvec[2*k+1] );
    psi.replaceGlobalValue( k, 0, z );
  }
}
// =============================================================================
// converts a real-valued vector to a complex-valued psi vector
void GlSystem::psi2real ( const Tpetra::MultiVector<double_complex,int> &psi,
                          Epetra_Vector                                 &x   )
{
  Teuchos::ArrayRCP<const double_complex> psiView = psi.getVector(0)->get1dView();
  for ( int k=0; k<NumComplexUnknowns; k++ )
    {
      x[2*k]   = real ( psiView[k] );
      x[2*k+1] = imag ( psiView[k] );
    }
}
// =============================================================================
void GlSystem::makeRealMap( const Teuchos::RCP<const Tpetra::Map<int> >  complexMap,
                            Teuchos::RCP<Epetra_Map>                     &realMap     )
{ 
  int numRealGlobalElements = 2*complexMap->getNodeNumElements();
  
  int myPid = TComm->getRank();

  // treat the phase condition on the first node
  if ( myPid==0 )
    numRealGlobalElements++;
  
  Epetra_IntSerialDenseVector   realMapGIDs( numRealGlobalElements );
  Teuchos::ArrayView<const int> myGlobalElements = complexMap->getNodeElementList();
  // Construct the map in such a way that all complex entries on processor K
  // are split up into real and imaginary part, which will both reside on
  // processor K again.
  for (unsigned int i=0; i<complexMap->getNodeNumElements(); i++)
    {
      realMapGIDs[2*i]   = 2*myGlobalElements[i];
      realMapGIDs[2*i+1] = 2*myGlobalElements[i] + 1;
    }
    
  // set the phase condition
  if ( myPid==0 )
    realMapGIDs[numRealGlobalElements-1] = NumRealUnknowns-1;

  realMap = Teuchos::rcp(new Epetra_Map( numRealGlobalElements,
                                         realMapGIDs.Length(),
                                         realMapGIDs.Values(),
                                         complexMap->getIndexBase(),
                                         *EComm                            ) );

  return;
}
// =============================================================================
bool GlSystem::computeF ( const Epetra_Vector &x,
                          Epetra_Vector       &FVec,
                          const NOX::Epetra::Interface::Required::FillType fillFlag )
{
//   vector<double_complex> psi ( NumComplexUnknowns );
//   double_complex         val;
  
//   // define the Tpetra platform
// #ifdef TPETRA_MPI
//         Tpetra::MpiPlatform<int, double_complex> platformV(MPI_COMM_WORLD);
//         Tpetra::MpiPlatform<int, int>            platformE(MPI_COMM_WORLD);
// #else
//         Tpetra::SerialPlatform<int, double_complex> platformV;
//         Tpetra::SerialPlatform<int, int>            platformE;
// #endif
  
  // ***************************************************************************
   
  // make sure that the input and output vectors are correctly mapped
  if ( !x.Map().SameAs( *RealMap ) ) {
      throw glException ( "GlSystem::computeF",
                          "Maps of x and the computed real-valued map do not coincide." );    
  }
  
  if ( !FVec.Map().SameAs( *RealMap ) ) {
      throw glException ( "GlSystem::computeF",
                          "Maps of FVec and the computed real-valued map do not coincide." );    
  }
   
  // define vector
  // @TODO: Replace this by Tpetra::Vector as soon as upstream is ready.
  Tpetra::MultiVector<double_complex,int> psi( ComplexMap, 1, true );

  // convert from x to psi2
  real2psi ( x, psi );

  // define output vector
  Tpetra::MultiVector<double_complex,int> res( ComplexMap, 1, true );
  
  // compute the GL residual
  res = Gl.computeGlVector ( psi );
  
  // transform back to fully real equation
  psi2real( res, FVec );
  
  // add phase condition
  FVec[2*NumComplexUnknowns] = 0.0;
  
//   // ***************************************************************************  
// 
//   // scatter x over all processors
//   Epetra_Export Exporter ( *StandardMap, *EverywhereMap );
//   xEverywhere.Export ( x, Exporter, Insert );
//   ( void ) real2complex ( xEverywhere, psi );
// 
//   // loop over the system rows
//   double passVal;
//   for ( int i=0; i<NumMyElements; i++ )
//     {
//       int myGlobalIndex = StandardMap->GID ( i );
//       if ( myGlobalIndex==2*NumComplexUnknowns )   // phase condition
//         {
//           passVal = 0.0;
//         }
//       else   // GL equations
//         {
//           // get the index of the complex valued equation
//           int psiIndex = realIndex2complexIndex ( myGlobalIndex );
//           // get the complex value
//           // TODO: The same value is actually fetched twice in this loop
//           //       possibly consecutively: Once for the real, once for
//           //       the imaginary part of it.
//           val = Gl.computeGl ( psiIndex, psi );
// 
// 
//           if ( ! ( myGlobalIndex%2 ) ) // myGlobalIndex is even
//             passVal = real ( val );
//           else // myGlobalIndex is odd
//             passVal = imag ( val );
//         }
//       FVec.ReplaceGlobalValues ( 1, &passVal, &myGlobalIndex );
//     }

  return true;
}
// =============================================================================
bool GlSystem::computeJacobian ( const Epetra_Vector &x,
                                 Epetra_Operator     &Jac )
{
  // compute the values of the Jacobian
  createJacobian ( VALUES, x );

  // optimize storage
  jacobian->FillComplete();

  // Sync up processors to be safe
  EComm->Barrier();

  return true;
}
// =============================================================================
bool GlSystem::computePreconditioner ( const Epetra_Vector    &x,
                                       Epetra_Operator        &Prec,
                                       Teuchos::ParameterList *precParams )
{
  throw glException ( "GlSystem::preconditionVector",
                      "Use explicit Jacobian only for this test problem!" );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector> GlSystem::getSolution()
{
  return initialSolution;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix> GlSystem::getJacobian()
{
  return jacobian;
}
/* =============================================================================
// Of an equation system
// \f[
// A\psi + B \psi^* = b
// \f]
// where \f$A,B\in\mathbb{C}^{n\times n}\f$, \f$\psi, b\in\mathbb{C}^{n}\f$,
// this routine constructs the corresponding real-valued equation system
// \f[
// \begin{pmatrix}
// \Re{A}+\Re{B} & -\Im{A}+\Im{B}\\
// \Im{A}+\Im{B} &  \Re{A}-\Re{B}
// \end{pmatrix}
// \begin{pmatrix}
// \Re{\psi}\\
// \Im{\psi}
// \end{pmatrix}
// =
// \begin{pmatrix}
// \Re{b}\\
// \Im{b}
// \end{pmatrix}
// \f].
// It also incorporates a phase condition.
*/
bool GlSystem::createJacobian ( const jacCreator    jc,
                                const Epetra_Vector &x )
{
  std::vector<double_complex> psi ( NumComplexUnknowns );

  vector<int>            colIndA, colIndB;
  vector<double_complex> valuesA, valuesB;

  int    *colInd      = NULL,
         *colIndAReal = NULL, *colIndAImag = NULL,
         *colIndBReal = NULL, *colIndBImag = NULL;
  double *values      = NULL,
         *valuesAReal = NULL, *valuesAImag = NULL,
         *valuesBReal = NULL, *valuesBImag = NULL;

  int ierr, k, complexRow, numEntries;
  
  if ( jc==VALUES )
    {
      Epetra_Vector xEverywhere ( *EverywhereMap );
      // scatter x over all processors
      Epetra_Export Exporter ( *RealMap, *EverywhereMap );

      xEverywhere.Export ( x, Exporter, Insert );
      real2complex ( xEverywhere, psi );
      // set the matrix to 0
      jacobian->PutScalar ( 0.0 );
    }
  else
    {
      if ( Graph.is_valid_ptr() )
        {
	  // Nullify Graph pointer.
          Graph = Teuchos::ENull();
        }
      // allocate the graph
      int approxNumEntriesPerRow = 1;
      Graph = Teuchos::rcp<Epetra_CrsGraph>( new Epetra_CrsGraph ( Copy, *RealMap, approxNumEntriesPerRow, false ) );
    }

  // Construct the Epetra Matrix
  for ( int i=0 ; i<NumMyElements ; i++ )
    {
      int Row = RealMap->GID ( i );

      if ( Row==2*NumComplexUnknowns )   // phase condition
        {
          // fill in phase condition stuff
          numEntries = 2*NumComplexUnknowns;
          colInd = new int[numEntries];
          for ( int k=0; k<numEntries; k++ )
            colInd[k] = k;
          if ( jc==VALUES )   // fill on columns and values
            {
              values = new double[numEntries];
              for ( int k=0; k<NumComplexUnknowns; k++ )
                {
                  values[2*k]   = -imag ( psi[k] );
                  values[2*k+1] =  real ( psi[k] );
                }
              // fill it in!
              ierr = jacobian->SumIntoGlobalValues ( Row,
                                                     numEntries,
                                                     values,
                                                     colInd );
              delete [] values;
              values = NULL;
            }
          else   // only fill the sparsity graph
            {
              Graph->InsertGlobalIndices ( Row,
                                           numEntries,
                                           colInd );
            }
          delete [] colInd;
          colInd = NULL;


        }
      else   // GL equations
        {
          // get the values and column indices
          // TODO: The same value is actually fetched twice in this loop
          //       possibly consecutively: Once for the real, once for
          //       the imaginary part of it.
          if ( ! ( Row%2 ) ) // Row even
            complexRow = Row/2;
          else
            complexRow = ( Row-1 ) /2;

          if ( jc==VALUES )
            {
              // fill on columns and values
              Gl.getJacobianRow ( complexRow,
                                  psi,
                                  colIndA, valuesA,
                                  colIndB, valuesB );
            }
          else
            {
              // only fill the sparsity graph
              Gl.getJacobianRowSparsity ( complexRow,
                                          colIndA,
                                          colIndB );
            }

          if ( ! ( Row%2 ) )   // myGlobalIndex is even <=> Real part of the equation system
            {
              // ---------------------------------------------------------------
              // insert the coefficients Re(A) of Re(psi)
              numEntries = colIndA.size();
              colIndAReal = new int[numEntries];
              for ( k=0; k<numEntries; k++ )
                colIndAReal[k] = 2*colIndA[k];

              if ( jc==VALUES )
                {
                  valuesAReal = new double[numEntries];
                  for ( k=0; k<numEntries; k++ )
                    valuesAReal[k] = real ( valuesA[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         numEntries,
                                                         valuesAReal,
                                                         colIndAReal );
                  delete [] valuesAReal;
                  valuesAReal = NULL;
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               numEntries,
                                               colIndAReal );
                }
              delete [] colIndAReal;
              colIndAReal = NULL;

              // insert the coefficients Re(B) of Re(psi)
              numEntries = colIndB.size();
              colIndBReal = new int[numEntries];
              for ( k=0; k<numEntries; k++ )
                colIndBReal[k] = 2*colIndB[k];
              if ( jc==VALUES )
                {
                  valuesBReal = new double[numEntries];
                  for ( k=0; k<numEntries; k++ )
                    valuesBReal[k] = real ( valuesB[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         numEntries,
                                                         valuesBReal,
                                                         colIndBReal );
                  delete [] valuesBReal;
                  valuesBReal = NULL;
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               numEntries,
                                               colIndBReal );
                }
              delete [] colIndBReal;
              colIndBReal = NULL;

              // insert the coefficients -Im(A) of Im(psi)
              numEntries = colIndA.size();
              colIndAImag = new int[numEntries];
              for ( k=0; k<numEntries; k++ )
                colIndAImag[k] = 2*colIndA[k]+1;
              if ( jc==VALUES )
                {
                  valuesAImag = new double[numEntries];
                  for ( k=0; k<numEntries; k++ )
                    valuesAImag[k] = -imag ( valuesA[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         numEntries,
                                                         valuesAImag,
                                                         colIndAImag );
                  delete [] valuesAImag;
                  valuesAImag = NULL;
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               numEntries,
                                               colIndAImag );
                }
              delete [] colIndAImag;
              colIndAImag = NULL;


              // insert the coefficients Im(B) of Im(psi)
              numEntries = colIndB.size();
              colIndBImag = new int[numEntries];
              for ( k=0; k<numEntries; k++ )
                colIndBImag[k] = 2*colIndB[k]+1;
              if ( jc==VALUES )
                {
                  valuesBImag = new double[numEntries];
                  for ( k=0; k<numEntries; k++ )
                    valuesBImag[k] = imag ( valuesB[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         numEntries,
                                                         valuesBImag,
                                                         colIndBImag );
                  delete [] valuesBImag;
                  valuesBImag = NULL;
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               numEntries,
                                               colIndBImag );
                }
              delete [] colIndBImag;
              colIndBImag = NULL;


              // right bordering
              int k = realIndex2complexIndex ( Row );
              int column = 2*NumComplexUnknowns;
              if ( jc==VALUES )
                {
                  double value  = imag ( psi[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         1,
                                                         &value,
                                                         &column );
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               1,
                                               &column );
                }

              // ---------------------------------------------------------
            }
          else   // Row is odd <=> Imaginary part of the equation
            {
              // ---------------------------------------------------------
              // insert the coefficients Im(A) of Re(psi)
              numEntries = colIndA.size();
              colIndAReal = new int[numEntries];
              for ( k=0; k<numEntries; k++ )
                colIndAReal[k] = 2*colIndA[k];
              if ( jc==VALUES )
                {
                  valuesAImag = new double[numEntries];
                  for ( k=0; k<numEntries; k++ )
                    valuesAImag[k] = imag ( valuesA[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         numEntries,
                                                         valuesAImag,
                                                         colIndAReal );
                  delete [] valuesAImag;
                  valuesAImag = NULL;
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               numEntries,
                                               colIndAReal );
                }
              delete [] colIndAReal;
              colIndAReal = NULL;

              // insert the coefficients Im(B) of Re(psi)
              numEntries = colIndB.size();
              colIndBReal = new int[numEntries];
              for ( k=0; k<numEntries; k++ )
                colIndBReal[k] = 2*colIndB[k];
              if ( jc==VALUES )
                {
                  valuesBImag = new double[numEntries];
                  for ( k=0; k<numEntries; k++ )
                    valuesBImag[k] = imag ( valuesB[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         numEntries,
                                                         valuesBImag,
                                                         colIndBReal );
                  delete [] valuesBImag;
                  valuesBImag = NULL;
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               numEntries,
                                               colIndBReal );
                }
              delete [] colIndBReal;
              colIndBReal = NULL;


              // insert the coefficients Re(A) of Im(psi)
              numEntries = colIndA.size();
              colIndAImag = new int[numEntries];
              for ( k=0; k<numEntries; k++ )
                colIndAImag[k] = 2*colIndA[k]+1;
              if ( jc==VALUES )
                {
                  valuesAReal = new double[numEntries];
                  for ( k=0; k<numEntries; k++ )
                    valuesAReal[k] = real ( valuesA[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         numEntries,
                                                         valuesAReal,
                                                         colIndAImag );
                  delete [] valuesAReal;
                  valuesAReal = NULL;
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               numEntries,
                                               colIndAImag );
                }
              delete [] colIndAImag;
              colIndAImag = NULL;


              // insert the coefficients -Re(B) of Im(psi)
              numEntries = colIndB.size();
              colIndBImag = new int[numEntries];
              for ( k=0; k<numEntries; k++ )
                colIndBImag[k] = 2*colIndB[k]+1;
              if ( jc==VALUES )
                {
                  valuesBReal = new double[numEntries];
                  for ( k=0; k<numEntries; k++ )
                    valuesBReal[k] = -real ( valuesB[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         numEntries,
                                                         valuesBReal,
                                                         colIndBImag );
                  delete [] valuesBReal;
                  valuesBReal = NULL;
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               numEntries,
                                               colIndBImag );
                }
              delete [] colIndBImag;
              colIndBImag = NULL;

              // right bordering
              int column = 2*NumComplexUnknowns;
              if ( jc==VALUES )
                {
                  int k = realIndex2complexIndex ( Row );
                  double value  = -real ( psi[k] );
                  ierr = jacobian->SumIntoGlobalValues ( Row,
                                                         1,
                                                         &value,
                                                         &column );
                }
              else
                {
                  Graph->InsertGlobalIndices ( Row,
                                               1,
                                               &column );
                }
              // ---------------------------------------------------------
            }
        }
    }

  // ---------------------------------------------------------------------------
  // finish up the graph construction
  try
    {
      if ( jc==VALUES )
        {
          jacobian->FillComplete();
          jacobian->OptimizeStorage();
        }
      else
        {
          Graph->FillComplete();
        }
    }
  catch ( int i )
    {
      std::string message = "FillComplete returned error code "
                            + EpetraExt::toString ( i ) + ".";
      throw glException ( "GlSystem::createJacobian",
                          message );
    }
  // ---------------------------------------------------------------------------

  // Sync up processors for safety's sake
  EComm->Barrier();
  
  return true;
}
// =============================================================================
// function used by LOCA
void GlSystem::setParameters ( const LOCA::ParameterVector &p )
{
  double h0 = p.getValue ( "H0" );

  // set H0 in the underlying problem class
  Gl.getStaggeredGrid()->setH0 ( h0 );
}
// =============================================================================
// function used by LOCA
void GlSystem::printSolution ( const Epetra_Vector &x,
                               double              conParam )
{
  static int conStep = 0;

  if (reverse)
      conStep--;
  else
      conStep++;

  // define vector
  // @TODO: Replace this by Tpetra::Vector as soon as upstream is ready.
  Tpetra::MultiVector<double_complex,int>  psi(ComplexMap,1,true);
  // convert from x to psi
  real2psi ( x, psi );

  // create temporary parameter list
  // TODO: get rid of this
  // -- An ugly thing here is that we have to explicitly mention the parameter
  // names. A solution could possibly be to include the parameter list in the
  // constructor of glSystem.
  Teuchos::ParameterList tmpList;
  tmpList.get ( "H0",         conParam );
  tmpList.get ( "edgelength", Gl.getStaggeredGrid()->getEdgelength() );
  tmpList.get ( "Nx",         Gl.getStaggeredGrid()->getNx() );
  tmpList.get ( "FE",         Gl.freeEnergy( psi ));

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // TODO:
  // see if the following can be replaced by some FORMAT construction
  std::string stepString;
  if (conStep>0)
      stepString = "+" + EpetraExt::toString ( conStep );
  else
      stepString = EpetraExt::toString ( conStep );
  std::string fileName = "continuationStep"
                         + stepString
                         + ".vtk";

  IoVirtual* fileIo = IoFactory::createFileIo ( outputDir+"/"+fileName );
  fileIo->write ( psi,
                  tmpList,
                  * ( Gl.getStaggeredGrid() ) );
  delete fileIo;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // fill the continuation parameters file
  std::string contFileBaseName;
  if (conStep>0)
       contFileBaseName = "continuation+.dat";
  else
       contFileBaseName = "continuation-.dat";
  std::string contFileName     = outputDir+"/"+contFileBaseName;

  std::ofstream contFileStream;

  // Set the output format
  // Think about replacing this with NOX::Utils::Sci.
  contFileStream.setf( std::ios::scientific );
  contFileStream.precision(15);

  if ( abs(conStep)==1 ) {
      contFileStream.open (contFileName.c_str(),ios::trunc);
      contFileStream << "# Step  "
                     << "\tH0              "
                     << "\tenergy              "
                     << "\t#vortices\n";
  } else {
      // just append to the the contents to the file
      contFileStream.open (contFileName.c_str(),ios::app);
  }

  contFileStream << "  " << conStep << "     "
                 << "\t" << conParam
                 << "\t" << Gl.freeEnergy   ( psi )
                 << "\t" << Gl.countVortices( psi ) << std::endl;

  contFileStream.close();
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}
// =============================================================================
// function used by LOCA
void GlSystem::setOutputDir ( const string &directory )
{
  outputDir = directory;
}
// =============================================================================
void GlSystem::solutionToFile ( const Epetra_Vector    &x,
                                Teuchos::ParameterList &problemParams,
                                const std::string      &fileName )
{
  
  // define vector
  // @TODO: Replace this by Tpetra::Vector as soon as upstream is ready.
  Tpetra::MultiVector<double_complex,int> psi(ComplexMap,1,true);

  // convert from x to psi2
  real2psi ( x, psi );

  problemParams.set( "FE", Gl.freeEnergy( psi ) );
  IoVirtual* fileIo = IoFactory::createFileIo ( fileName );
  fileIo->write ( psi,
                  problemParams,
                  * ( Gl.getStaggeredGrid() ) );

}
// =============================================================================