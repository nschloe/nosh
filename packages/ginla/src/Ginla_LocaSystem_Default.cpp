#include "Ginla_LocaSystem_Default.h"

#include "Ginla_Helpers.h"
#include "Ginla_Operator_Virtual.h"
#include "Ginla_Perturbation_Virtual.h"
#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_StatsWriter.h"

#include <Epetra_Map.h>
#include <LOCA_Stepper.H>
#include <Teuchos_DefaultComm.hpp>
#include <Epetra_CrsMatrix.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// =============================================================================
// Default constructor without perturbation
Ginla::LocaSystem::Default::
Default ( const Teuchos::RCP<Ginla::Operator::Virtual>           & glOperator,
          const Teuchos::RCP<const Epetra_Comm>                  & eComm,
          const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > & complexMap,
          const Teuchos::RCP<Ginla::IO::StatsWriter>             & statsWriter,
          const Teuchos::RCP<Ginla::IO::StateWriter>             & stateWriter
        ) :
        glOperator_ ( glOperator ),
        perturbation_ ( Teuchos::null ),
        preconditioner_( Teuchos::null ),
        stepper_ ( Teuchos::null ),
        komplex_ ( Teuchos::rcp ( new Ginla::Komplex ( eComm, complexMap ) ) ),
        statsWriter_( statsWriter ),
        stateWriter_( stateWriter ),
        firstTime_ ( true )
{
}
// =============================================================================
// Destructor
Ginla::LocaSystem::Default::~Default()
{
    stepper_ = Teuchos::null;
    return;
}
// =============================================================================
Teuchos::RCP<Ginla::State>
Ginla::LocaSystem::Default::
createState( const Epetra_Vector & x
           ) const
{
    const Teuchos::RCP<ComplexVector> psi = komplex_->real2complex ( x );
    return Teuchos::rcp( new Ginla::State( psi, glOperator_->getGrid() ) );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::LocaSystem::Default::
createX(  const Ginla::State & state )
{
    return komplex_->complex2real ( state.getValuesConst() );
}
// =============================================================================
bool
Ginla::LocaSystem::Default::
computeF ( const Epetra_Vector & x,
           Epetra_Vector       & FVec,
           const NOX::Epetra::Interface::Required::FillType fillFlag )
{ 
    // make sure that the input and output vectors are correctly mapped
    TEST_FOR_EXCEPTION ( !x.Map().SameAs ( *komplex_->getRealMap() ),
                         std::logic_error,
                         "Maps of x and the computed real-valued map do not coincide." );

    TEST_FOR_EXCEPTION ( !FVec.Map().SameAs ( *komplex_->getRealMap() ),
                         std::logic_error,
                         "Maps of FVec and the computed real-valued map do not coincide." );
    // ------------------------------------------------------------------------
    // convert from x to a state
    const Teuchos::RCP<const Ginla::State> state = this->createState( x );
    // ------------------------------------------------------------------------
    // compute the GL residual

    Teuchos::RCP<Ginla::State> res = glOperator_->getF( state );

    // add perturbation
    if ( !perturbation_.is_null() )
    {
        Teuchos::ArrayRCP<double_complex> resView = res->getValuesNonConst()->get1dViewNonConst();
        // loop over the nodes
        for ( unsigned int k=0; k<state->getValuesConst()->getLocalLength(); k++ )
        {
            int globalIndex = state->getValuesConst()->getMap()->getGlobalElement ( k );
            resView[k] += perturbation_->computePerturbation ( globalIndex );
        }
    }
    // ------------------------------------------------------------------------
    // TODO Avoid this explicit copy.
    // transform back to fully real equation
    FVec = *(this->createX( *res ));
    // ------------------------------------------------------------------------

    return true;
}
// =============================================================================
bool
Ginla::LocaSystem::Default::
computeJacobian ( const Epetra_Vector & x,
                  Epetra_Operator     & Jac )
{
    Teuchos::Array<int> indicesA, indicesB;
    Teuchos::Array<double_complex> valuesA, valuesB;

    komplex_->zeroOutMatrix();

    const Teuchos::RCP<const Ginla::State> state = this->createState( x );
//     Teuchos::RCP<ComplexVector> psi = komplex_->real2complex ( x );

    // update to the latest psi vector before retrieving the Jacobian
//     glOperator_->updatePsi ( psi );
    
    // loop over the rows and fill the matrix
    int numMyElements = komplex_->getComplexMap()->getNodeNumElements();
    for ( int row = 0; row < numMyElements; row++ )
    {
        int globalRow = komplex_->getComplexMap()->getGlobalElement ( row );
        // get the values from the operator
        glOperator_->getJacobianRow ( state,
                                      globalRow,
                                      indicesA, valuesA,
                                      indicesB, valuesB );
        // ... and fill them into glKomplex_
        komplex_->updateRow ( globalRow,
                              indicesA, valuesA,
                              indicesB, valuesB,
                              firstTime_ );
    }

    if ( firstTime_ )
    {
        komplex_->finalizeMatrix();
        firstTime_ = false;
    }

    return true;
}
// =============================================================================
bool
Ginla::LocaSystem::Default::
computePreconditioner ( const Epetra_Vector & x,
                        Epetra_Operator     & Prec,
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
Teuchos::RCP<Epetra_CrsMatrix>
Ginla::LocaSystem::Default::getJacobian() const
{
    return komplex_->getMatrix();
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
Ginla::LocaSystem::Default::getPreconditioner() const
{
    return preconditioner_;
}
// =============================================================================
bool
Ginla::LocaSystem::Default::
computeShiftedMatrix ( double alpha,
                       double beta,
                       const Epetra_Vector & x,
                       Epetra_Operator     & A )
{
    // compute the values of the Jacobian
    computeJacobian ( x, A );

    komplex_->getMatrix()->Scale ( alpha );

    Epetra_Vector newDiag ( x );
    Epetra_Vector unitVector ( x );
    unitVector.PutScalar ( 1.0 );
    //  newDiag.PutScalar(0.0);
    komplex_->getMatrix()->ExtractDiagonalCopy ( newDiag );
    newDiag.Update ( beta, unitVector, 1.0 );
    komplex_->getMatrix()->ReplaceDiagonalValues ( newDiag );

    return true;
}
// =============================================================================
// function used by LOCA
void
Ginla::LocaSystem::Default::
setParameters ( const LOCA::ParameterVector & p )
{
    TEUCHOS_ASSERT( glOperator_.is_valid_ptr() && !glOperator_.is_null() );
    glOperator_->setParameters ( p );
    
    if ( perturbation_.is_valid_ptr() && !perturbation_.is_null() )
        perturbation_->setParameters ( p );

    return;
}
// =============================================================================
void
Ginla::LocaSystem::Default::
setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper )
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
Ginla::LocaSystem::Default::
releaseLocaStepper()
{
    stepper_ = Teuchos::null;
}
// =============================================================================
// function used by LOCA
void
Ginla::LocaSystem::Default::
printSolution ( const  Epetra_Vector &x,
                double conParam
              )
{
    // define vector
    const Teuchos::RCP<const Ginla::State> state = this->createState( x );

    // The switch hack is necessary as different continuation algorithms
    // call printSolution() a different number of times per step, e.g.,
    // to store solutions, null vectors, and so forth.
    switch ( continuationType_ )
    {
    case ONEPARAMETER:
        this->printSolutionOneParameterContinuation ( state );
        break;
    case TURNINGPOINT:
        this->printSolutionTurningPointContinuation ( state );
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
printSolution ( const Epetra_Vector & x,
                const std::string   & filenameAppendix
              ) const
{
    const Teuchos::RCP<const Ginla::State> state = this->createState( x );

    stateWriter_->write( state,
                         stepper_->getStepNumber(),
                         "-" + filenameAppendix
                       );
    return;
}
// =============================================================================
void
Ginla::LocaSystem::Default::
printSolutionOneParameterContinuation ( const Teuchos::RCP<const Ginla::State> & state )
{
    stateWriter_->write( state,
                         stepper_->getStepNumber(),
                         *(glOperator_->getParameters())
                         );

    writeContinuationStats ( state );
    return;
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
printSolutionTurningPointContinuation ( const Teuchos::RCP<const Ginla::State> & state )
{
    static bool isSolution = false;

    // alternate between solution and nullvector
    isSolution = !isSolution;

    // determine file name
    stringstream baseName;
    if ( isSolution )
    {
        stateWriter_->write( state,
                             stepper_->getStepNumber(),
                             "-state",
                             *(glOperator_->getParameters())
                           );
        writeContinuationStats ( state );
    }
    else
    {
        stateWriter_->write( state,
                             stepper_->getStepNumber(),
                             "-nullvector"
                           );
    }

    return;
}
// =============================================================================
void
Ginla::LocaSystem::Default::
writeContinuationStats ( const Teuchos::RCP<const Ginla::State> & state )
{   
    TEUCHOS_ASSERT( statsWriter_.is_valid_ptr()
                    && !statsWriter_.is_null() );

    TEUCHOS_ASSERT( statsWriter_->getList().is_valid_ptr()
                    && !statsWriter_->getList().is_null() );

    Teuchos::RCP<Teuchos::ParameterList> paramList = statsWriter_->getList();
                    
    paramList->set( "0step", stepper_->getStepNumber() );
    paramList->set( "2#nonlinear steps", stepper_->getSolver()->getNumIterations() );
  
    // put the parameter list into statsWriter_
    std::string labelPrepend = "1";
    Ginla::Helpers::appendToTeuchosParameterList( *(statsWriter_->getList()),
                                                  *(glOperator_->getParameters()),
                                                  labelPrepend );

    TEUCHOS_ASSERT( state.is_valid_ptr() && !state.is_null() );
    
    paramList->set( "2free energy", state->freeEnergy() );
    paramList->set( "2||x||_2 scaled", state->normalizedScaledL2Norm() );
    paramList->set( "2vorticity", state->getVorticity() );

    // actually print the data
    statsWriter_->print();

    return;
}
// =============================================================================
// function used by LOCA
void
Ginla::LocaSystem::Default::
setOutputDir ( const string & directory )
{
    stateWriter_->setOutputDir( directory );
    return;
}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
Ginla::LocaSystem::Default::
create_CommInt ( const Teuchos::RCP<const Epetra_Comm> & epetraComm )
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
        Teuchos::RCP<const Teuchos::SerialComm<int> >
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
Ginla::LocaSystem::Default::
getMap() const
{
    return komplex_->getRealMap();
}
// =============================================================================
Teuchos::RCP<const Ginla::Komplex>
Ginla::LocaSystem::Default::
getKomplex() const
{
    return komplex_;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::LocaSystem::Default::
createSystemVector( const Teuchos::ParameterList & p )
{
   TEST_FOR_EXCEPTION( true,
                       std::logic_error,
                       "Not yet implemented." );
}
// =============================================================================