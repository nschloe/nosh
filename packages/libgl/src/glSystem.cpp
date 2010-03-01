#include "glSystem.h"

#include <Epetra_Export.h>
#include <Epetra_CrsMatrix.h>
#include <NOX_Utils.H>

#include <EpetraExt_RowMatrixOut.h>

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
GlSystem::GlSystem ( GinzburgLandau::GinzburgLandau &gl,
                     const Teuchos::RCP<const Epetra_Comm> eComm,
                     const Teuchos::RCP<const ComplexVector> psi,
                     const std::string outputDir,
                     const std::string outputDataFileName,
                     const std::string outputFileFormat,
                     const std::string solutionFileNameBase,
                     const std::string nullvectorFileNameBase,
                     const unsigned int maxNumDigits
                   ) :
        stepper_ ( Teuchos::null ),
        glKomplex_ ( Teuchos::rcp ( new GlKomplex ( eComm,psi->getMap() ) ) ),
        Gl_ ( gl ),
        preconditioner_( Teuchos::null ),
        initialSolution_ ( glKomplex_->complex2real ( psi ) ),
        outputDir_ ( outputDir ),
        solutionFileNameBase_ ( solutionFileNameBase ),
        nullvectorFileNameBase_ ( nullvectorFileNameBase ),
        outputFileFormat_ ( outputFileFormat ),
        outputDataFileName_ ( outputDataFileName ),
        firstTime_ ( true ),
        maxNumDigits_( maxNumDigits )
{
}
// =============================================================================
// Destructor
GlSystem::~GlSystem()
{
    stepper_ = Teuchos::null;
}
// =============================================================================
bool
GlSystem::computeF ( const Epetra_Vector &x,
                     Epetra_Vector &FVec,
                     const NOX::Epetra::Interface::Required::FillType fillFlag )
{
    // make sure that the input and output vectors are correctly mapped
    TEST_FOR_EXCEPTION ( !x.Map().SameAs ( *glKomplex_->getRealMap() ),
                         std::logic_error,
                         "Maps of x and the computed real-valued map do not coincide." );

    TEST_FOR_EXCEPTION ( !FVec.Map().SameAs ( *glKomplex_->getRealMap() ),
                         std::logic_error,
                         "Maps of FVec and the computed real-valued map do not coincide." );
                         
    // convert from x to psi
    const Teuchos::RCP<ComplexVector> psi = glKomplex_->real2complex ( x );

    // compute the GL residual
    Teuchos::RCP<ComplexVector> res = Gl_.computeGlVector ( psi );

    // TODO Avoid this explicit copy.
    // transform back to fully real equation
    FVec = * ( glKomplex_->complex2real ( *res ) );

    return true;
}
// =============================================================================
bool
GlSystem::computeJacobian ( const Epetra_Vector &x, Epetra_Operator &Jac )
{
    // compute the values of the Jacobian
    createJacobian ( x );

    return true;
}
// =============================================================================
bool GlSystem::computePreconditioner ( const Epetra_Vector &x,
                                       Epetra_Operator &Prec,
                                       Teuchos::ParameterList *precParams )
{
//  Epetra_Vector diag = x;
//  diag.PutScalar(1.0);
//  preconditioner_->ReplaceDiagonalValues( diag );

    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Use explicit Jacobian only for this test problem!" );
    return true;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
GlSystem::getSolution() const
{
    return initialSolution_;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GlSystem::getJacobian() const
{
//      TEUCHOS_ASSERT( jacobian_.is_valid_ptr() && !jacobian_.is_null() );
    return glKomplex_->getMatrix();
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GlSystem::getPreconditioner() const
{
    return preconditioner_;
}
// =============================================================================
void
GlSystem::createJacobian ( const Epetra_Vector &x )
{
    Teuchos::Array<int> indicesA, indicesB;
    Teuchos::Array<double_complex> valuesA, valuesB;

    glKomplex_->zeroOutMatrix();

    Teuchos::RCP<ComplexVector> psi = glKomplex_->real2complex ( x );

    // loop over the rows and fill the matrix
    int numMyElements = glKomplex_->getComplexMap()->getNodeNumElements();
    for ( int row = 0; row < numMyElements; row++ )
    {
        int globalRow = glKomplex_->getComplexMap()->getGlobalElement ( row );
        // get the values from Gl_
        Gl_.getJacobianRow ( globalRow, psi, indicesA, valuesA, indicesB, valuesB );
        // ... and fill them into glKomplex_
        glKomplex_->updateRow ( globalRow, indicesA, valuesA, indicesB, valuesB, firstTime_ );
    }

    if ( firstTime_ )
    {
        glKomplex_->finalizeMatrix();
        firstTime_ = false;
    }

    return;
}
// =============================================================================
bool
GlSystem::computeShiftedMatrix ( double alpha,
                                 double beta,
                                 const Epetra_Vector &x,
                                 Epetra_Operator &A )
{
    // compute the values of the Jacobian
    createJacobian ( x );

    glKomplex_->getMatrix()->Scale ( alpha );

    Epetra_Vector newDiag ( x );
    Epetra_Vector unitVector ( x );
    unitVector.PutScalar ( 1.0 );
    //  newDiag.PutScalar(0.0);
    glKomplex_->getMatrix()->ExtractDiagonalCopy ( newDiag );
    newDiag.Update ( beta, unitVector, 1.0 );
    glKomplex_->getMatrix()->ReplaceDiagonalValues ( newDiag );

    return true;
}
// =============================================================================
// function used by LOCA
void
GlSystem::setParameters ( const LOCA::ParameterVector & p )
{

    TEST_FOR_EXCEPTION ( !p.isParameter ( "H0" ),
                         std::logic_error,
                         "Label \"H0\" not valid." );
    double h0 = p.getValue ( "H0" );
    Gl_.setH0 ( h0 );

    TEST_FOR_EXCEPTION ( !p.isParameter ( "scaling" ),
                         std::logic_error,
                         "Label \"scaling\" not valid." );
    double scaling = p.getValue ( "scaling" );
    Gl_.setScaling ( scaling );

    if ( p.isParameter ( "chi" ) )
    {
        double chi = p.getValue ( "chi" );
        Gl_.setChi ( chi );
    }
}
// =============================================================================
void
GlSystem::setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper )
{
    stepper_ = stepper;

    // extract the continuation type
    const Teuchos::ParameterList & bifurcationSublist =
        stepper_->getList()->sublist ( "LOCA" ).sublist ( "Bifurcation" );

    std::string bifurcationType = bifurcationSublist.get<string> ( "Type" );

    if ( bifurcationType == "None" )
        continuationType_ = ONEPARAMETER;
    else if ( bifurcationType == "Turning Point" )
        continuationType_ = TURNINGPOINT;
    else
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Unknown continuation type \""
                             << bifurcationType << "\"." );
}
// =============================================================================
void
GlSystem::releaseLocaStepper()
{
    stepper_ = Teuchos::null;
}
// =============================================================================
// function used by LOCA
void
GlSystem::printSolution ( const Epetra_Vector &x,
                          double conParam )
{
    // define vector
    const Teuchos::RCP<ComplexVector> psi =
        glKomplex_->real2complex ( x );

    // The switch hack is necessary as different continuation algorithms
    // call printSolution() a different number of times per step, e.g.,
    // to store solutions, null vectors, and so forth.
    switch ( continuationType_ )
    {
    case ONEPARAMETER:
        printSolutionOneParameterContinuation ( psi );
        break;
    case TURNINGPOINT:
        printSolutionTurningPointContinuation ( psi );
        break;
    default:
        TEST_FOR_EXCEPTION ( true,
                             std::logic_error,
                             "Illegal continuation type " << continuationType_ );
    }
}
// =============================================================================
void
GlSystem::printSolutionOneParameterContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                                ) const
{
    static int conStep = -1;
    conStep++;
    
    stringstream fileName;
    fileName
    << outputDir_ << "/" << solutionFileNameBase_ 
    << setw ( maxNumDigits_ ) << setfill ( '0' ) << conStep << ".vtk";
    
    // actually print the state to fileName
    Gl_.writeSolutionToFile ( psi, fileName.str() );

    writeContinuationStats ( conStep, psi );
}
// =============================================================================
// In Turning Point continuation, the printSolution method is called exactly
// twice per step:
//
//   1. For printing the solution.
//   2. For printing the right null vector of the Jacobian.
//
// The method gets called subsequently in this order.
void
GlSystem::printSolutionTurningPointContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                                ) const
{
    static bool printSolution=false;
    static int conStep = -1;

    // alternate between solution and nullvector
    printSolution = !printSolution;

    // increment the step counter only when printing a solution
    if ( printSolution )
        conStep++;

    // determine file name
    stringstream fileName;
    if ( printSolution )
    {
        fileName
        << outputDir_ << "/" << solutionFileNameBase_
        << setw ( maxNumDigits_ ) << setfill ( '0' ) << conStep << ".vtk";
        writeContinuationStats ( conStep, psi );
    }
    else
        fileName
        << outputDir_ << "/" << nullvectorFileNameBase_
        << setw ( maxNumDigits_ ) << setfill ( '0' ) << conStep << ".vtk";

    // actually print the state to fileName
    Gl_.writeSolutionToFile ( psi, fileName.str() );
}
// =============================================================================
void
GlSystem::writeContinuationStats ( const int conStep,
                                   const Teuchos::RCP<const ComplexVector> psi ) const
{
    // fill the continuation parameters file
    std::string contFileName = outputDir_ + "/" + outputDataFileName_;
    std::ofstream contFileStream;

    // Set the output format
    // Think about replacing this with NOX::Utils::Sci.
    contFileStream.setf ( std::ios::scientific );
    contFileStream.precision ( 15 );

    if ( conStep == 0 )
    {
        contFileStream.open ( contFileName.c_str(), ios::trunc );
        contFileStream << "# Step  \t";
        Gl_.appendStats ( contFileStream, true );
        contFileStream << "\t#nonlinear steps\n";
    }
    else
    {
        // just append to the the contents to the file
        contFileStream.open ( contFileName.c_str(), ios::app );
    }

    int nonlinearIterations = stepper_->getSolver()->getNumIterations();

    contFileStream << "  " << conStep << "      \t";
    Gl_.appendStats ( contFileStream, false, psi );
    contFileStream << "       \t" << nonlinearIterations << std::endl;

    contFileStream.close();
}
// =============================================================================
// function used by LOCA
void GlSystem::setOutputDir ( const string &directory )
{
    outputDir_ = directory;
}
// =============================================================================
void
GlSystem::writeSolutionToFile ( const Epetra_Vector & x,
                                const std::string   & filePath
                              ) const
{
    // TODO: Remove the need for several real2complex calls per step.
    Gl_.writeSolutionToFile ( glKomplex_->real2complex ( x ), filePath );
}
// =============================================================================
void
GlSystem::writeAbstractStateToFile ( const Epetra_Vector & x,
                                     const std::string   & filePath
                                   ) const
{
    // TODO: Remove the need for several real2complex calls per step.
    Gl_.writeAbstractStateToFile ( glKomplex_->real2complex ( x ), filePath );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
GlSystem::getGlSystemVector ( const Teuchos::RCP<const ComplexVector> psi ) const
{
    return glKomplex_->complex2real ( *psi );
}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
GlSystem::create_CommInt ( const Teuchos::RCP<const Epetra_Comm> & epetraComm )
{
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::set_extra_data;

#ifdef HAVE_MPI
    RCP<const Epetra_MpiComm>
    mpiEpetraComm = rcp_dynamic_cast<const Epetra_MpiComm> ( epetraComm );
    if ( mpiEpetraComm.get() )
    {
        RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
        rawMpiComm = Teuchos::opaqueWrapper ( mpiEpetraComm->Comm() );
        set_extra_data ( mpiEpetraComm, "mpiEpetraComm", Teuchos::inOutArg ( rawMpiComm ) );
        RCP<const Teuchos::MpiComm<int> >
        mpiComm = rcp ( new Teuchos::MpiComm<int> ( rawMpiComm ) );
        return mpiComm;
    }
#else
    RCP<const Epetra_SerialComm>
    serialEpetraComm = rcp_dynamic_cast<const Epetra_SerialComm> ( epetraComm );
    if ( serialEpetraComm.get() )
    {
        RCP<const Teuchos::SerialComm<int> >
        serialComm = rcp ( new Teuchos::SerialComm<int>() );
        set_extra_data ( serialEpetraComm, "serialEpetraComm", Teuchos::inOutArg ( serialComm ) );
        return serialComm;
    }
#endif // HAVE_MPI

    // If you get here then the conversion failed!
    return Teuchos::null;
}
// =============================================================================
double
GlSystem::getH0() const
{
    return Gl_.getH0();
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
GlSystem::getRealMap() const
{
    return glKomplex_->getRealMap();
}
// =============================================================================
Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> >
GlSystem::getComplexMap() const
{
    return glKomplex_->getComplexMap();
}
// =============================================================================
void
GlSystem::setH0 ( const double h0 )
{
    Gl_.setH0 ( h0 );
}
// =============================================================================
void
GlSystem::setScaling ( const double scaling )
{
    Gl_.setScaling ( scaling );
}
// =============================================================================
void
GlSystem::setChi ( const double chi )
{
    Gl_.setChi ( chi );
}
// =============================================================================
Teuchos::RCP<const GlKomplex>
GlSystem::getGlKomplex() const
{
    return glKomplex_;
}
// =============================================================================
