#include "Ginla_LocaSystem_Default.h"
#include "Ginla_Helpers.h"

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
Ginla::LocaSystem::Default::
Default ( GinzburgLandau::GinzburgLandau          & gl,
          const Teuchos::RCP<const Epetra_Comm>     eComm,
          const Teuchos::RCP<const ComplexVector>   psi,
          const std::string outputDir,
          const std::string outputDataFileName,
          const std::string solutionFileNameBase,
          const unsigned int maxNumDigits
        ) :
        stepper_ ( Teuchos::null ),
        glKomplex_ ( Teuchos::rcp ( new Ginla::Komplex ( eComm,psi->getMap() ) ) ),
        Gl_ ( gl ),
        preconditioner_( Teuchos::null ),
        initialSolution_ ( glKomplex_->complex2real ( psi ) ),
        outputDir_ ( outputDir ),
        solutionFileNameBase_ ( solutionFileNameBase ),
        outputDataFileName_ ( outputDataFileName ),
        firstTime_ ( true ),
        maxNumDigits_( maxNumDigits )
{
  return;
}
// =============================================================================
// Destructor
Ginla::LocaSystem::Default::~Default()
{
    stepper_ = Teuchos::null;
//     contFileStream.close();
  return;
}
// =============================================================================
bool
Ginla::LocaSystem::Default::
computeF ( const Epetra_Vector &x,
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
Ginla::LocaSystem::Default::
computeJacobian ( const Epetra_Vector & x,
                  Epetra_Operator     & Jac )
{
    // compute the values of the Jacobian
    createJacobian ( x );

    return true;
}
// =============================================================================
bool Ginla::LocaSystem::Default::computePreconditioner ( const Epetra_Vector &x,
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
Ginla::LocaSystem::Default::getSolution() const
{
    return initialSolution_;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
Ginla::LocaSystem::Default::getJacobian() const
{
//      TEUCHOS_ASSERT( jacobian_.is_valid_ptr() && !jacobian_.is_null() );
    return glKomplex_->getMatrix();
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
Ginla::LocaSystem::Default::getPreconditioner() const
{
    return preconditioner_;
}
// =============================================================================
void
Ginla::LocaSystem::Default::createJacobian ( const Epetra_Vector &x )
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
Ginla::LocaSystem::Default::
computeShiftedMatrix ( double alpha,
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
Ginla::LocaSystem::Default::setParameters ( const LOCA::ParameterVector & p )
{
  Gl_.setParameters( p );
}
// =============================================================================
void
Ginla::LocaSystem::Default::setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper )
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
Ginla::LocaSystem::Default::releaseLocaStepper()
{
    stepper_ = Teuchos::null;
}
// =============================================================================
// function used by LOCA
void
Ginla::LocaSystem::Default::printSolution ( const Epetra_Vector &x,
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
Ginla::LocaSystem::Default::
printSolutionOneParameterContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                      )
{
    int conStep = stepper_->getStepNumber();
    
    stringstream baseName;
    baseName
    << outputDir_ << "/" << solutionFileNameBase_ 
    << setw ( maxNumDigits_ ) << setfill ( '0' ) << conStep;
        
    // actually print the state to fileName
    Ginla::Helpers::writeStateToFile( psi,
                                      Gl_.getOperator()->getGrid(),
                                      *(Gl_.getOperator()->getParameters()),
                                      baseName.str(),
                                      "VTI" );

    writeContinuationStats ( psi );
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
Ginla::LocaSystem::Default::
printSolutionTurningPointContinuation ( const Teuchos::RCP<const ComplexVector> & psi
                                      )
{
    static bool printSolution=false;
    int conStep = stepper_->getStepNumber();

    // alternate between solution and nullvector
    printSolution = !printSolution;

    // determine file name
    stringstream baseName;
    if ( printSolution )
    {
        baseName
        << outputDir_ << "/" << solutionFileNameBase_
        << setw ( maxNumDigits_ ) << setfill ( '0' ) << conStep
        << "-state";
        writeContinuationStats ( psi );
    }
    else
        baseName
        << outputDir_ << "/" << solutionFileNameBase_
        << setw ( maxNumDigits_ ) << setfill ( '0' ) << conStep
        << "-nullvector";
    
    // actually print the state to fileName
    Ginla::Helpers::writeStateToFile( psi,
                                      Gl_.getOperator()->getGrid(),
                                      *(Gl_.getOperator()->getParameters()),
                                      baseName.str(),
                                      "VTI" );

    return;
}
// =============================================================================
void
Ginla::LocaSystem::Default::
writeContinuationStats ( const Teuchos::RCP<const ComplexVector> & psi )
{   
    TEUCHOS_ASSERT( Gl_.getStatsWriter().is_valid_ptr()
                    && !Gl_.getStatsWriter().is_null() );

    TEUCHOS_ASSERT( Gl_.getStatsWriter()->getList().is_valid_ptr()
                    && !Gl_.getStatsWriter()->getList().is_null() );

    Teuchos::RCP<Teuchos::ParameterList> paramList = Gl_.getStatsWriter()->getList();
                    
    paramList->set( "0step", stepper_->getStepNumber() );
    paramList->set( "2#nonlinear steps", stepper_->getSolver()->getNumIterations() );
  
    // put the parameter list into statsWriter_
    std::string labelPrepend = "1";
    Ginla::Helpers::appendToTeuchosParameterList( *(Gl_.getStatsWriter()->getList()),
                                                  *(Gl_.getOperator()->getParameters()),
                                                  labelPrepend );

    TEUCHOS_ASSERT( psi.is_valid_ptr() && !psi.is_null() );
    
    paramList->set( "2free energy", Ginla::Helpers::freeEnergy ( *psi, *(Gl_.getOperator()->getGrid()) ) );
    paramList->set( "2||x||_2 scaled", Ginla::Helpers::normalizedScaledL2Norm ( *psi, *(Gl_.getOperator()->getGrid()) ) );
    paramList->set( "2vorticity", Ginla::Helpers::getVorticity ( *psi, *(Gl_.getOperator()->getGrid()) ) );

    // actually print the data
    Gl_.getStatsWriter()->print();

    return;
}
// =============================================================================
// function used by LOCA
void Ginla::LocaSystem::Default::setOutputDir ( const string &directory )
{
    outputDir_ = directory;
}
// =============================================================================
void
Ginla::LocaSystem::Default::
writeSolutionToFile ( const Epetra_Vector & x,
                      const std::string   & filePath
                    ) const
{   
    // TODO: Remove the need for several real2complex calls per step.
    Ginla::Helpers::writeStateToFile( glKomplex_->real2complex ( x ),
                                      Gl_.getOperator()->getGrid(),
                                      *(Gl_.getOperator()->getParameters()),
                                      filePath,
                                      "VTI" );
    return;
}
// =============================================================================
void
Ginla::LocaSystem::Default::
writeAbstractStateToFile ( const Epetra_Vector & x,
                           const std::string   & filePath
                         ) const
{
    // TODO: Remove the need for several real2complex calls per step.
    LOCA::ParameterVector empty;
    Ginla::Helpers::writeStateToFile( glKomplex_->real2complex ( x ),
                                      Gl_.getOperator()->getGrid(),
                                      empty,
                                      filePath,
                                      "VTI" );
    return;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::LocaSystem::Default::getGlSystemVector ( const Teuchos::RCP<const ComplexVector> psi ) const
{
    return glKomplex_->complex2real ( *psi );
}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
Ginla::LocaSystem::Default::create_CommInt ( const Teuchos::RCP<const Epetra_Comm> & epetraComm )
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
Teuchos::RCP<const Epetra_Map>
Ginla::LocaSystem::Default::getRealMap() const
{
    return glKomplex_->getRealMap();
}
// =============================================================================
Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> >
Ginla::LocaSystem::Default::getComplexMap() const
{
    return glKomplex_->getComplexMap();
}
// =============================================================================
Teuchos::RCP<const Ginla::Komplex>
Ginla::LocaSystem::Default::getGlKomplex() const
{
    return glKomplex_;
}
// =============================================================================
