#include "glPredictorSystem.h"
#include "ioFactory.h"

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

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// abbreviate the complex type name
typedef std::complex<double> double_complex;

// =============================================================================
// Default constructor
GlPredictorSystem::GlPredictorSystem( GinzburgLandau::GinzburgLandau &gl,
	                                  const Teuchos::RCP<const Epetra_Comm> eComm,
	                                  const Teuchos::RCP<ComplexVector> psi,
	                                  const Teuchos::RCP<ComplexVector> tangent,
	                                  const Teuchos::RCP<ComplexVector> predictor,
	                                  const std::string outputDir,
                                      const std::string outputDataFileName,
                                      const std::string outputFileFormat,
	                                  const std::string solutionFileNameBase,
	                                  const std::string nullvectorFileNameBase
                                    ) :
	glSystem_( Teuchos::rcp(new GlSystemWithConstraint( gl,
                                          eComm,
                                          psi,
                                          outputDir,
                                          outputDataFileName,
                                          outputFileFormat,
                                          solutionFileNameBase,
                                          nullvectorFileNameBase
                                         ) ) ),
   EComm_( eComm ),
   psi_ ( psi ),
   tangent_ ( tangent ),
   predictor_ ( predictor ),
   preconditioner_( 0 )
{
  // TODO Don't throw exception in constructor?
  TEST_FOR_EXCEPTION( !tangent.is_valid_ptr(),
                      std::logic_error,
	                  "Invalid pointer" );

  // TODO Don't throw exception in constructor?
  TEST_FOR_EXCEPTION( tangent.is_null(),
                      std::logic_error,
	                  "Input guess is null pointer" );

  sliceMap_    = glSystem_->getMap();
  extendedMap_ = createExtendedMap( *sliceMap_ );
  cout << "GGG DONE THE CONSTRUCTOR" << endl;

}
// =============================================================================
// Destructor
GlPredictorSystem::~GlPredictorSystem() {
}
// =============================================================================
bool
GlPredictorSystem::computeF( const Epetra_Vector &xExtended,
		                     Epetra_Vector &FVecExtended,
		                     const NOX::Epetra::Interface::Required::FillType fillFlag)
{
	cout << "11" << endl;

  // make sure that the input and output vectors are correctly mapped
  TEST_FOR_EXCEPTION( !xExtended.Map().SameAs(*extendedMap_),
                      std::logic_error,
	                  "Maps of x and the computed real-valued map do not coincide." );

  TEST_FOR_EXCEPTION( !FVecExtended.Map().SameAs(*extendedMap_),
                      std::logic_error,
                      "Maps of FVec and the computed real-valued map do not coincide." );

  cout << "AAA" << endl;

  Epetra_Vector x( *sliceMap_ );
  for (int k=0; k<x.MyLength(); k++ ) {
    int kGlobal = extendedMap_->GID(k);
    x.ReplaceMyValue( k,0, xExtended[kGlobal] );
  }

  cout << "bbb" << endl;

  Epetra_Vector FVec( *sliceMap_ );
  glSystem_->computeF( x, FVec );

  cout << "CCC" << endl;

  for (int k=0; k < FVec.MyLength(); k++ ) {
    int kGlobal = FVec.Map().GID(k);
    FVecExtended(kGlobal) = FVec(kGlobal);
  }

  cout << "DDD" << endl;

  Teuchos::RCP<const Epetra_Vector> tangentReal = glSystem_->getGlKomplex()->complex2real(*tangent_);
  Teuchos::RCP<Epetra_Vector> predictorReal     = glSystem_->getGlKomplex()->complex2real(*predictor_);
  Teuchos::RCP<Epetra_Vector> psiReal           = glSystem_->getGlKomplex()->complex2real(*psi_);

  Teuchos::RCP<Epetra_Vector> diff = Teuchos::rcp( new Epetra_Vector( psiReal->Map() ) );

  cout << "EEE" << endl;

  // compute last entry
  for (int k=0; k<diff->MyLength(); k++ )
	  diff->ReplaceMyValue( k, 0, (*predictorReal)[k] - (*psiReal)[k] );

  int k = FVecExtended.GlobalLength()-1;
  double *result;
  tangentReal->Dot( *diff,result );
  FVecExtended[k] = result[0];

  cout << "FFF" << endl;

  return true;
}
// =============================================================================
bool GlPredictorSystem::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac) {

	cout << "GGG" << endl;

  // compute the values of the Jacobian
  createJacobian(x);

  // optimize storage
  jacobian_->FillComplete();

  // Sync up processors to be safe
  EComm_->Barrier();

  return true;
}
// =============================================================================
bool GlPredictorSystem::computePreconditioner( const Epetra_Vector &x,
	                                                 Epetra_Operator &Prec,
                                                     Teuchos::ParameterList *precParams )
{
	cout << "22" << endl;
    TEST_FOR_EXCEPTION( true,
			std::logic_error,
	                "Use explicit Jacobian only for this test problem!" );
  return true;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
GlPredictorSystem::getSolution() const {
	cout << "33" << endl;
	return initialSolution_;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GlPredictorSystem::getJacobian() const {
	cout << "44" << endl;
	cout << jacobian_.is_valid_ptr() << endl;
	cout << jacobian_.is_null() << endl;
	cout << "44b" << endl;
	return jacobian_;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GlPredictorSystem::getPreconditioner() const {
	cout << "55" << endl;
	cout << preconditioner_.is_valid_ptr() << endl;
	cout << preconditioner_.is_null() << endl;
    cout << "55b" << endl;
    return preconditioner_;
}
// =============================================================================
bool GlPredictorSystem::createJacobian(const Epetra_Vector &xExtended)
{

	cout << "6" << endl;
		// approximate df/dH0 by finite difference
		double origH0 = glSystem_->getH0();

		Epetra_Vector x( *sliceMap_ );
		  for (int k=0; k<x.MyLength(); k++ ) {
		    int kGlobal = extendedMap_->GID(k);
		    x.ReplaceMyValue( k,0, xExtended[kGlobal] );
		  }

		Epetra_Vector FVec( *sliceMap_ );
		Epetra_Vector FVec2( *sliceMap_ );
		glSystem_->computeF( x, FVec );
        // reset parameters
		double eps = 1.0e-6;
		glSystem_->setH0( origH0+eps );
		glSystem_->computeF( x, FVec2 );

		Epetra_Vector finiteDiff( *sliceMap_ );
		for(int k=0; k<finiteDiff.MyLength(); k++ ) {
			int kGlobal = finiteDiff.Map().GID(k);
			finiteDiff[kGlobal] = (FVec2[kGlobal] - FVec[kGlobal])
					            / eps;
		}
		// set back to original value
		glSystem_->setH0( origH0 );

		Teuchos::RCP<Epetra_CrsMatrix> originalJac = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *sliceMap_, 1, false ) );
		glSystem_->computeJacobian( x, *originalJac );

		originalJac = glSystem_->getJacobian();

		// loop through the matrix and insert all values originalJac into jacobian_
		int     numRowEntries;
		double* rowValues;
		int*    rowIndices;
		for(int row=0; row<x.GlobalLength(); row++) {
			originalJac->ExtractMyRowView( row, numRowEntries, rowValues, rowIndices );
			PutMyValues( jacobian_, row, numRowEntries, rowValues, rowIndices );

			int rowGlobal = finiteDiff.Map().GID(row);
			double* val;
			val[0] = finiteDiff[rowGlobal];
			int* col;
			col[0] = x.GlobalLength();
			PutMyValues( jacobian_, row, 1, val, col );
		}

		// add last row
		Teuchos::RCP<Epetra_Vector> xTangent = glSystem_->getGlKomplex()->complex2real(*tangent_);
		int lastRow = x.GlobalLength();
        for( int column=0; column<xTangent->GlobalLength(); column++) {
        	double* val;
        	val[0] = (*xTangent)[column];
        	int* col;
        	col[0] = column;
		    PutMyValues( jacobian_, lastRow, 1, val, col );
        }
        // don't explicitly fill 0.0 at (end,end)

	// ---------------------------------------------------------------------------
	// finish up the graph construction
	try {
		if (firstTime_) {
			jacobian_->FillComplete();
			jacobian_->OptimizeStorage();
			firstTime_ = false;
		}
	} catch (int i) {
	    TEST_FOR_EXCEPTION( true,
				            std::logic_error,
				            "FillComplete returned error code " << i );
	}
	// ---------------------------------------------------------------------------

	// Sync up processors for safety's sake
	EComm_->Barrier();

	return true;
}
// =============================================================================
void
GlPredictorSystem::PutMyValues( Teuchos::RCP<Epetra_CrsMatrix> mat,
		                        int row,
		                        int numEntries,
		                        double* values,
		                        int* columns ) const
{
	cout << "77" << endl;
  if ( firstTime_ ) {
	  TEST_FOR_EXCEPT( mat->InsertMyValues ( row, numEntries, values, columns )<0 );
  } else {
	  TEST_FOR_EXCEPT( mat->ReplaceMyValues( row, numEntries, values, columns )<0 );
  }
}
// =============================================================================
void
GlPredictorSystem::PutGlobalValues( Teuchos::RCP<Epetra_CrsMatrix> mat,
		                            int     row,
		                            int     numEntries,
		                            double* values,
		                            int*    columns ) const
{
	cout << "8" << endl;
  if ( firstTime_ ) {
	  TEST_FOR_EXCEPT( mat->InsertGlobalValues ( row, numEntries, values, columns )<0 );
  } else {
	  TEST_FOR_EXCEPT( mat->ReplaceGlobalValues( row, numEntries, values, columns )<0 );
  }
}
// =============================================================================
// function used by LOCA
void
GlPredictorSystem::setParameters(const LOCA::ParameterVector &p) {

	cout << "HHH" << endl;

  TEST_FOR_EXCEPTION( !p.isParameter("H0"),
                      std::logic_error,
                      "Label \"H0\" not valid." );
  double h0 = p.getValue("H0");
  glSystem_->setH0(h0);

  TEST_FOR_EXCEPTION( !p.isParameter("scaling"),
                      std::logic_error,
                      "Label \"scaling\" not valid." );
  double scaling = p.getValue("scaling");
  glSystem_->setScaling( scaling );

  if (p.isParameter("chi")) {
      double chi = p.getValue("chi");
      glSystem_->setChi( chi );
  }

  cout << "III" << endl;
}
// =============================================================================
void
GlPredictorSystem::setLocaStepper(const Teuchos::RCP<const LOCA::Stepper> stepper)
{
	cout << "GGGLOCASTEPPER" << endl;
	glSystem_->setLocaStepper( stepper );
}
// =============================================================================
void
GlPredictorSystem::releaseLocaStepper()
{
	glSystem_->releaseLocaStepper();
}
// =============================================================================
// function used by LOCA
void
GlPredictorSystem::printSolution( const Epetra_Vector &xExtended,
                                  double conParam )
{
	cout << "KKK" << endl;

	Epetra_Vector xSliced( *sliceMap_ );
	for (int k=0; k<xSliced.MyLength(); k++ ) {
	  int kGlobal = extendedMap_->GID(k);
	  xSliced.ReplaceMyValue( k,0, xExtended[kGlobal] );
	}

	glSystem_->printSolution( xSliced, conParam );

//	// define vector
//	const Teuchos::RCP<ComplexVector> psi =
//			              Teuchos::rcp( new ComplexVector(ComplexMap_, true) );
//	// convert from x to psi
//	real2complex(x, *psi);
//
//	static int conStep = -1;
//	conStep++;
//
//	std::string fileName = outputDir_ + "/" + solutionFileNameBase_
//			+ EpetraExt::toString(conStep) + ".vtk";
//
//	// actually print the state to fileName
//	Gl_.writeSolutionToFile(psi, fileName);
//
//	writeContinuationStats( conStep, psi );
}
// =============================================================================
// function used by LOCA
void GlPredictorSystem::setOutputDir(const string &directory) {
	glSystem_->setOutputDir( directory );
	cout << "HIHI" << endl;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
GlPredictorSystem::createExtendedMap( const Epetra_BlockMap & sliceMap  )
{
  // fill up realMapGIDs
  int numMyElements = sliceMap.NumMyElements();
  Teuchos::Array<int> myElements( numMyElements );
  sliceMap.MyGlobalElements( myElements.getRawPtr() );

  // Unconditionally put the phase constraint on the first process.
  int myPID = sliceMap.Comm().MyPID();
  if ( myPID==0 ) {
    int n = sliceMap.NumGlobalElements();
    // extend the GIDs by the phase constraint
    myElements.append( n );
  }

  int numGlobalElements = sliceMap.NumGlobalElements() + 1;
  Teuchos::RCP<Epetra_Map> extendedMap = Teuchos::rcp( new Epetra_Map(numGlobalElements,
                                             myElements.length(),
                                             myElements.getRawPtr(),
                                             sliceMap.IndexBase(),
                                  			 sliceMap.Comm() )
                            );
  return extendedMap;
}
// =============================================================================
