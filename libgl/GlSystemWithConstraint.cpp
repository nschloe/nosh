/*
 * GlSystemWithConstraint.cpp
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schloemer
 */

#include "GlSystemWithConstraint.h"

#include "ioFactory.h"

#include <boost/format.hpp>

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
GlSystemWithConstraint::GlSystemWithConstraint( GinzburgLandau::GinzburgLandau &gl,
                                                const Teuchos::RCP<const Epetra_Comm> eComm,
                                                const Teuchos::RCP<ComplexVector> psi,
                                                const std::string outputDir,
                                                const std::string outputDataFileName,
                                                const std::string outputFileFormat,
                                                const std::string solutionFileNameBase,
                                                const std::string nullvectorFileNameBase,
                                                const unsigned int maxStepNumberDecimals
                                              ) :
        glSystem_( gl, eComm, psi, outputDir, outputDataFileName, outputFileFormat, solutionFileNameBase, nullvectorFileNameBase ),
        NumMyElements_(0),
        NumComplexUnknowns_(0),
        Gl_(gl),
        EComm_(eComm),
        regularRealMap_(0),
        extendedRealMap_(0),
        ComplexMap_(0),
        rhs_(0),
        Graph_(0),
        jacobian_(0),
        initialSolution_(0),
        solutionFileNameBase_(solutionFileNameBase),
        nullvectorFileNameBase_(nullvectorFileNameBase),
        outputFileFormat_(outputFileFormat),
        outputDataFileName_(outputDataFileName),
        glKomplex_( Teuchos::rcp(new GlKomplex(eComm,psi->getMap()) ) ),
        maxStepNumberDecimals_( maxStepNumberDecimals ),
        firstTime_(true)
{
  NumComplexUnknowns_ = Gl_.getNumUnknowns();

  // TODO Don't throw exception in constructor?
  TEST_FOR_EXCEPTION( psi->getGlobalLength() != (unsigned int) NumComplexUnknowns_,
                      std::logic_error,
                      "Size of the initial guess vector ("
                      << psi->getGlobalLength()
                      << ") does not coincide with the number of unknowns ("
                      << NumComplexUnknowns_ << ")" );

  ComplexMap_ = glSystem_.getComplexMap();

  // do the rest of the initialization
  initialize( psi );
}
// =============================================================================
// constructor *without* initial guess
GlSystemWithConstraint::GlSystemWithConstraint(GinzburgLandau::GinzburgLandau &gl,
                   const Teuchos::RCP<const Epetra_Comm> eComm,
                   const std::string outputDir,
                   const std::string outputDataFileName,
                   const std::string outputFileFormat,
                   const std::string solutionFileNameBase,
                   const std::string nullvectorFileNameBase,
                   const unsigned int maxStepNumberDecimals
                  ) :
glSystem_( gl, eComm, outputDir, outputDataFileName, outputFileFormat, solutionFileNameBase, nullvectorFileNameBase ),
NumMyElements_(0),
NumComplexUnknowns_(0),
Gl_(gl),
EComm_(eComm),
regularRealMap_(0),
extendedRealMap_(0),
ComplexMap_(0),
rhs_(0),
Graph_(0),
jacobian_(0),
initialSolution_(0),
solutionFileNameBase_(solutionFileNameBase),
nullvectorFileNameBase_(nullvectorFileNameBase),
outputFileFormat_(outputFileFormat),
outputDataFileName_(outputDataFileName),
glKomplex_( Teuchos::null ),
maxStepNumberDecimals_( maxStepNumberDecimals ),
firstTime_(true)
{
  NumComplexUnknowns_ = Gl_.getNumUnknowns();

  // TODO There is (until now?) no way to convert a Teuchos::Comm (of psi)
  // to an Epetra_Comm (of the real valued representation of psi), so the
  // Epetra_Comm has to be generated explicitly, and two communicators are kept
  // side by side all the time. One must make sure that the two are actually
  // equivalent, which can be checked by Thyra's conversion method create_Comm.
  // TODO Is is actually necessary to have equivalent communicators on the
  // real-valued and the complex-valued side?
  // How to compare two communicators anyway?

  // define complex map
  ComplexMap_ = glSystem_.getComplexMap();

  glKomplex_ = Teuchos::rcp(new GlKomplex(eComm,ComplexMap_) );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // initialize solution
  Teuchos::RCP<ComplexVector> psi = Teuchos::rcp(new ComplexVector(ComplexMap_));
  // TODO Move default initialization out to main file
  double_complex alpha(1.0, 0.0);
  psi->putScalar(alpha); // default initialization
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // do the rest of the initialization
  initialize( psi );
}
// =============================================================================
// Destructor
GlSystemWithConstraint::~GlSystemWithConstraint() {
}
// =============================================================================
void
GlSystemWithConstraint::initialize(const Teuchos::RCP<ComplexVector> psi)
{
	Teuchos::RCP<Epetra_Vector> tmp = glKomplex_->complex2real(*psi);

	// Create the maps with and without phase constraint.
	regularRealMap_  = glSystem_.getRealMap();
	extendedRealMap_ = createExtendedRealMap( *regularRealMap_ );

	initialSolution_ = Teuchos::rcp( new Epetra_Vector(*extendedRealMap_), true );
	for (int k=0; k<tmp->MyLength(); k++ ) {
	    initialSolution_->ReplaceMyValue( k, 0, (*tmp)[tmp->Map().GID(k)] );
	}
	int n = initialSolution_->GlobalLength();
	initialSolution_->ReplaceGlobalValue( n-1, 0, 0.0 );

	NumMyElements_ = extendedRealMap_->NumMyElements();

	// TODO Why is *this necessary??
	jacobian_ = Teuchos::rcp( new Epetra_CrsMatrix(Copy, *extendedRealMap_, 0) );

	// Initialize the format for the the continuation step number.
	// Here: 00012 for step no. 12, if maxStepNumberDecimals_=5.
    stepNumFileNameFormat_ = boost::str(  boost::format("%%|0%d|") % maxStepNumberDecimals_ );
}
// =============================================================================
bool
GlSystemWithConstraint::computeF(const Epetra_Vector &x,
                                       Epetra_Vector &FVec,
                                 const NOX::Epetra::Interface::Required::FillType fillFlag)
{
  TEST_FOR_EXCEPTION( !regularRealMap_.is_valid_ptr() || regularRealMap_.is_null(),
                      std::logic_error,
                      "regularRealMap_ not properly initialized." );

  TEST_FOR_EXCEPTION( !extendedRealMap_.is_valid_ptr() || extendedRealMap_.is_null(),
                      std::logic_error,
                      "extendedRealMap_ not properly initialized." );

  // make sure that the input and output vectors are correctly mapped
  TEST_FOR_EXCEPTION( !x.Map().SameAs(*extendedRealMap_),
                      std::logic_error,
                      "Maps of x and the computed real-valued map do not coincide. "
                      << "Check, for example, the number of elements "
                      << "(" << x.Map().NumGlobalElements() << " for x vs. "
                      << extendedRealMap_->NumGlobalElements() << " for extendedRealMap_).");

  TEST_FOR_EXCEPTION( !FVec.Map().SameAs(*extendedRealMap_),
                      std::logic_error,
                      "Maps of FVec and the computed real-valued map do not coincide."
                      << "Check, for example, the number of elements "
                      << "(" << FVec.Map().NumGlobalElements() << " for FVec vs. "
                      << extendedRealMap_->NumGlobalElements() << " for extendedRealMap_).");

  // TODO replace by {im,ex}porter
  // strip off the phase constraint
  Epetra_Vector tmp(*regularRealMap_);
  for (int k=0; k<tmp.MyLength(); k++)
    tmp.ReplaceMyValue( k, 0, x[x.Map().GID(k)] );

  Epetra_Vector shortFVec(*regularRealMap_);
  glSystem_.computeF( tmp, shortFVec, fillFlag );

  // copy over and add phase condition
  // TODO replace by {im,ex}porter
  for (int k=0; k<shortFVec.MyLength(); k++) {
    FVec.ReplaceMyValue( k, 0, shortFVec[shortFVec.Map().GID(k)] );
  }
  FVec.ReplaceGlobalValue( 2*NumComplexUnknowns_, 0, 0.0 );

  return true;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
GlSystemWithConstraint::createExtendedRealMap( const Epetra_BlockMap & realMap  ) const
{
  // fill up realMapGIDs
  int numMyElements = realMap.NumMyElements();
  Teuchos::Array<int> myElements( numMyElements );
  realMap.MyGlobalElements( myElements.getRawPtr() );

  // Unconditionally put the phase constraint on the first process.
  int myPID = realMap.Comm().MyPID();
  if ( myPID==0 ) {
    int n = realMap.NumGlobalElements();
    // extend the GIDs by the phase constraint
    myElements.append( n );
  }

  int numGlobalElements = realMap.NumGlobalElements() + 1;
  return Teuchos::rcp( new Epetra_Map(numGlobalElements,
                                      myElements.length(),
                                      myElements.getRawPtr(),
                                      realMap.IndexBase(),
                                      realMap.Comm() )
                     );
}
// =============================================================================
bool GlSystemWithConstraint::computeJacobian(const Epetra_Vector   & x,
		                                           Epetra_Operator & Jac
		                                    )
{
  // strip off the phase constraint
  Epetra_Vector tmp(*regularRealMap_);
  for (int k=0; k<tmp.MyLength(); k++)
    tmp.ReplaceMyValue( k, 0, x[x.Map().GID(k)] );

  // TODO Strip down Jac, too.
  // --   Not really necessary as it's not being used anyway.

  // compute the underlying Jacobian
  glSystem_.computeJacobian( tmp, Jac );

  // compute the values of the Jacobian
  createJacobian( x );

  // optimize storage
//  jacobian_->FillComplete();

  // Sync up processors to be safe
  EComm_->Barrier();

  return true;
}
// =============================================================================
bool GlSystemWithConstraint::computePreconditioner( const Epetra_Vector    & x,
                                                    Epetra_Operator        & Prec,
                                                    Teuchos::ParameterList * precParams )
{
  TEST_FOR_EXCEPTION( true,
                      std::logic_error,
                      "Use explicit Jacobian only for this test problem!" );
  return true;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
GlSystemWithConstraint::getSolution() const {
        return initialSolution_;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GlSystemWithConstraint::getJacobian() const {
        return jacobian_;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GlSystemWithConstraint::getPreconditioner() const {
        return preconditioner_;
}
// =============================================================================
// It also incorporates a phase condition.
void
GlSystemWithConstraint::createJacobian( const Epetra_Vector & x   )
{
  TEST_FOR_EXCEPTION( !extendedRealMap_.is_valid_ptr() || extendedRealMap_.is_null(),
                      std::logic_error,
                      "extendedRealMap_ not properly initialized." );

  // TODO replace by {im,ex}porter
  // strip off the phase constraint
  Epetra_Vector tmp(*regularRealMap_);
  for (int k=0; k<tmp.MyLength(); k++)
    tmp.ReplaceMyValue( k, 0, x[x.Map().GID(k)] );

  // TODO don't explicitly construct psi? get1dCopy on the rhs
  Teuchos::RCP<ComplexVector>             psi     = glKomplex_->real2complex(tmp);
  Teuchos::ArrayRCP<const double_complex> psiView = psi->get1dView();

  // get the unbordered Jacobian
  Teuchos::RCP<Epetra_CrsMatrix> regularJacobian = glSystem_.getJacobian();

  std::string::size_type n = regularRealMap_->NumGlobalElements();
  // create right bordering
  Teuchos::Array<double> rightBorder(n);
  for (int k = 0; k < NumComplexUnknowns_; k++) {
      rightBorder[2 * k]     =  imag(psiView[k]);
      rightBorder[2 * k + 1] = -real(psiView[k]);
  }

  // create lower border
  Teuchos::Array<double> lowerBorder(n);
  for (int k = 0; k < NumComplexUnknowns_; k++) {
   	  lowerBorder[2 * k]     = -imag(psiView[k]);
      lowerBorder[2 * k + 1] =  real(psiView[k]);
  }

  // corner element
  double d = 0.0;

  // create the bordered Jacobian out of this
  fillBorderedMatrix( jacobian_, regularJacobian, rightBorder, lowerBorder, d, firstTime_ );

  if (firstTime_) {
  	  TEST_FOR_EXCEPT( 0 != jacobian_->FillComplete() );
  	  TEST_FOR_EXCEPT( 0 != jacobian_->OptimizeStorage() );
  	  firstTime_ = false;
  }

  // Sync up processors for safety's sake
  EComm_->Barrier();

  return;
}
// =============================================================================
bool
GlSystemWithConstraint::computeShiftedMatrix( double alpha,
		                                      double beta,
                                              const Epetra_Vector &x,
                                              Epetra_Operator &A)
{
        // compute the values of the Jacobian
        createJacobian(x);

        jacobian_->Scale(alpha);
        //  jacobian_->FillComplete();

        Epetra_Vector newDiag(x);
        Epetra_Vector unitVector(x);
        unitVector.PutScalar(1.0);
        //  newDiag.PutScalar(0.0);
        jacobian_->ExtractDiagonalCopy(newDiag);
        newDiag.Update(beta, unitVector, 1.0);
        jacobian_->ReplaceDiagonalValues(newDiag);

        // Sync up processors to be safe
        EComm_->Barrier();

        return true;
}
// =============================================================================
// function used by LOCA
void
GlSystemWithConstraint::setParameters(const LOCA::ParameterVector &p) {

  TEST_FOR_EXCEPTION( !p.isParameter("H0"),
                      std::logic_error,
                      "Label \"H0\" not valid." );
  double h0 = p.getValue("H0");
  Gl_.setH0(h0);

  TEST_FOR_EXCEPTION( !p.isParameter("scaling"),
                      std::logic_error,
                      "Label \"scaling\" not valid." );
  double scaling = p.getValue("scaling");
  Gl_.setScaling( scaling );

  if (p.isParameter("chi")) {
      double chi = p.getValue("chi");
      Gl_.setChi( chi );
  }
}
// =============================================================================
void
GlSystemWithConstraint::setLocaStepper(const Teuchos::RCP<const LOCA::Stepper> stepper)
{
	glSystem_.setLocaStepper( stepper );
}
// =============================================================================
void
GlSystemWithConstraint::releaseLocaStepper()
{
	glSystem_.releaseLocaStepper();
}
// =============================================================================
// function used by LOCA
void
GlSystemWithConstraint::printSolution( const  Epetra_Vector &x,
                                       double conParam )
{
	glSystem_.printSolution(  x, conParam );
}
// =============================================================================
// function used by LOCA
void GlSystemWithConstraint::setOutputDir(const string &directory) {
    glSystem_.setOutputDir( directory );
}
// =============================================================================
void
GlSystemWithConstraint::writeSolutionToFile( const Epetra_Vector &x,
                                             const std::string &filePath) const
{
	glSystem_.writeSolutionToFile( x, filePath );
}
// =============================================================================
void
GlSystemWithConstraint::writeAbstractStateToFile( const Epetra_Vector &x,
                                                  const std::string &filePath) const
{
	glSystem_.writeAbstractStateToFile( x, filePath );
}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
GlSystemWithConstraint::create_CommInt( const Teuchos::RCP<const Epetra_Comm> &epetraComm )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

#ifdef HAVE_MPI
  RCP<const Epetra_MpiComm>
    mpiEpetraComm = rcp_dynamic_cast<const Epetra_MpiComm>(epetraComm);
  if( mpiEpetraComm.get() ) {
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = Teuchos::opaqueWrapper(mpiEpetraComm->Comm());
    set_extra_data( mpiEpetraComm, "mpiEpetraComm", Teuchos::inOutArg(rawMpiComm) );
    RCP<const Teuchos::MpiComm<int> >
      mpiComm = rcp(new Teuchos::MpiComm<int>(rawMpiComm));
    return mpiComm;
  }
#else
  RCP<const Epetra_SerialComm>
    serialEpetraComm = rcp_dynamic_cast<const Epetra_SerialComm>(epetraComm);
  if( serialEpetraComm.get() ) {
    RCP<const Teuchos::SerialComm<int> >
      serialComm = rcp(new Teuchos::SerialComm<int>());
    set_extra_data( serialEpetraComm, "serialEpetraComm", Teuchos::inOutArg(serialComm) );
    return serialComm;
  }
#endif // HAVE_MPI

  // If you get here then the conversion failed!
  return Teuchos::null;
}
// =============================================================================
// TODO delete?
const Teuchos::RCP<const GlKomplex>
GlSystemWithConstraint::getGlKomplex() const {
  return glKomplex_;
}
// =============================================================================
double
GlSystemWithConstraint::getH0() const
{
  return Gl_.getH0();
}
// =============================================================================
const Teuchos::RCP<const Epetra_Map>
GlSystemWithConstraint::getMap() const
{
  return extendedRealMap_;
}
// =============================================================================
void
GlSystemWithConstraint::setH0(const double h0)
{
  Gl_.setH0( h0 );
}
// =============================================================================
void
GlSystemWithConstraint::setScaling(const double scaling)
{
  Gl_.setScaling( scaling );
}
// =============================================================================
void
GlSystemWithConstraint::setChi(const double chi)
{
  Gl_.setChi( chi );
}
// =============================================================================
void
GlSystemWithConstraint::fillBorderedMatrix( const Teuchos::RCP<      Epetra_CrsMatrix> & extendedMatrix,
                                            const Teuchos::RCP<const Epetra_CrsMatrix> & regularMatrix,
                                            const Teuchos::Array<double>               & rightBorder,
                                            // TODO Declare the following const as soon as Trilinos allows (ReplaceGlobalValues)
                                                  Teuchos::Array<double>               & lowerBorder,
                                                  double                                 d,
                                                  bool                                   firstTime
                                          ) const
{
  TEUCHOS_ASSERT( regularMatrix.is_valid_ptr()  && !regularMatrix.is_null() );
  TEUCHOS_ASSERT( extendedMatrix.is_valid_ptr() && !extendedMatrix.is_null() );

  int m = regularMatrix->NumGlobalRows();
  int n = regularMatrix->NumGlobalCols();

  // check if the sizes all match
  TEUCHOS_ASSERT_EQUALITY( m+1, extendedMatrix->NumGlobalRows() );
  TEUCHOS_ASSERT_EQUALITY( n+1, extendedMatrix->NumGlobalCols() );
  TEUCHOS_ASSERT_EQUALITY( m, rightBorder.length() );
  TEUCHOS_ASSERT_EQUALITY( n, lowerBorder.length() );

  int numMyRows = regularMatrix->NumMyRows();

  // fill the matrix with the entries
  int numRowNonZeros;

  int maxNumEntries = regularMatrix->MaxNumEntries() + 1; // count the last column in
  int    * indices = new int   [maxNumEntries];
  double * values  = new double[maxNumEntries];
  for ( int myRow=0; myRow<numMyRows; myRow++ ) {
	  // extract row view
	  TEUCHOS_ASSERT_EQUALITY( 0, regularMatrix->ExtractMyRowView( myRow, numRowNonZeros, values, indices ) );

	  // Can't use InsertMyIndices because the *indices are given in global indexing.
	  int globalRow = extendedMatrix->Map().GID(myRow);

      // write the data to the new matrix
	  TEUCHOS_ASSERT_INEQUALITY( 0, <=, PutRow( extendedMatrix, globalRow, numRowNonZeros, values, indices, firstTime ) );

      // add last column
      double val = rightBorder[globalRow];
      TEUCHOS_ASSERT_INEQUALITY( 0, <=, PutRow( extendedMatrix, globalRow, 1, &val, &n, firstTime ) );
  }

  // set the last row
  // create the indices array
  int lastRow[n];
  for ( int k=0; k<n; k++ )
	  lastRow[k] = k;
  TEUCHOS_ASSERT_INEQUALITY( 0, <=, PutRow( extendedMatrix, n, n, lowerBorder.getRawPtr(), lastRow, firstTime ) );

  // set the last element d
  TEUCHOS_ASSERT_INEQUALITY( 0, <=, PutRow( extendedMatrix, n, 1, &d, &n, firstTime ) );

  return;
}
// =============================================================================
int
GlSystemWithConstraint::PutRow( const Teuchos::RCP<Epetra_CrsMatrix> A,
		                        int      Row,
              		            const int      numIndices,
		                        double * values,
		                        int    * indices,
		                        const bool     firstTime ) const
{
	if (firstTime) {
		return A->InsertGlobalValues( Row, numIndices, values, indices );
	} else {
		return A->ReplaceGlobalValues( Row, numIndices, values, indices );
	}
}
// =============================================================================
