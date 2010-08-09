#include "Ginla_FDM_LocaSystem_Default.h"

#include "Ginla_Helpers.h"
#include "Ginla_FDM_Operator_Virtual.h"
#include "Ginla_FDM_Perturbation_Virtual.h"
#include "Ginla_IO_StateWriter.h"
#include "Ginla_IO_StatsWriter.h"
#include "Ginla_Komplex_DoubleMatrix.h"

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
Ginla::FDM::LocaSystem::Default::
Default ( const Teuchos::RCP<Ginla::FDM::Operator::Virtual>      & glOperator,
          const Teuchos::RCP<const Epetra_Comm>                  & eComm,
          const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > & complexMap,
          const Teuchos::RCP<Ginla::IO::StatsWriter>             & statsWriter,
          const Teuchos::RCP<Ginla::IO::StateWriter>             & stateWriter
        ) :
        glOperator_ ( glOperator ),
        perturbation_ ( Teuchos::null ),
        preconditioner_( Teuchos::null ),
        stepper_ ( Teuchos::null ),
        komplex_ ( Teuchos::rcp ( new Ginla::Komplex::LinearProblem ( eComm, complexMap ) ) ),
        statsWriter_( statsWriter ),
        stateWriter_( stateWriter ),
        firstTime_ ( true )
{
}
// =============================================================================
// Destructor
Ginla::FDM::LocaSystem::Default::~Default()
{
    stepper_ = Teuchos::null;
    return;
}
// =============================================================================
bool
Ginla::FDM::LocaSystem::Default::
computeF ( const Epetra_Vector & x,
           Epetra_Vector       & FVec,
           const NOX::Epetra::Interface::Required::FillType fillFlag
         )
{ 
    // make sure that the input and output vectors are correctly mapped
    TEST_FOR_EXCEPTION ( !x.Map().SameAs ( *komplex_->getRealMap() ),
                         std::logic_error,
                         "Maps of x and the computed real-valued map do not coincide." );

    TEST_FOR_EXCEPTION ( !FVec.Map().SameAs ( *komplex_->getRealMap() ),
                         std::logic_error,
                         "Maps of FVec and the computed real-valued map do not coincide."
                       );
    // ------------------------------------------------------------------------
    // convert from x to a state
    const Teuchos::RCP<const Ginla::FDM::State> state =
        this->createFdmState_( x );
    // ------------------------------------------------------------------------
    // compute the GL residual
    Teuchos::RCP<Ginla::State::Virtual> res = glOperator_->getF( state );

    // add perturbation
    if ( !perturbation_.is_null() )
    {
        Teuchos::ArrayRCP<double_complex> resView =
            res->getPsiNonConst()->get1dViewNonConst();
        // loop over the nodes
        for ( unsigned int k=0; k<state->getPsi()->getLocalLength(); k++ )
        {
            int globalIndex = state->getPsi()->getMap()->getGlobalElement ( k );
            resView[k] += perturbation_->computePerturbation ( globalIndex );
        }
    }
    // ------------------------------------------------------------------------
    // TODO Avoid this explicit copy.
    // transform back to fully real equation
    FVec = *( this->createSystemVector( *res ) );
    // ------------------------------------------------------------------------
    
    return true;
}
// =============================================================================
bool
Ginla::FDM::LocaSystem::Default::
computeJacobian ( const Epetra_Vector & x,
                  Epetra_Operator     & Jac )
{ 
    // In Ginla::Komplex, the values are merely sumIntoLocalValue'd,
    // so make sure we set this to zero from the start.
    komplex_->zeroOutMatrix();

    const Teuchos::RCP<const Ginla::FDM::State> state =
        this->createFdmState_( x );
    
    // create a real-valued matrix of the AB-Jacobian
    komplex_->update( glOperator_->getJacobian( state ), firstTime_ );
    
    Teuchos::RCP<const Ginla::Komplex::DoubleMatrix> AB = glOperator_->getJacobian( state );
    Teuchos::RCP<const ComplexMatrix> A = AB->getMatrixA();
    Teuchos::RCP<const ComplexMatrix> B = AB->getMatrixB();
    
//   std::cout << "AAA" << std::endl;
//   Teuchos::ArrayRCP<const long int> indices;
//   Teuchos::ArrayRCP<const double_complex> values;
//   for ( long int i=0; i<A->getRowMap()->getGlobalNumElements(); i++ )
//   {
//       A->getLocalRowView( i, indices, values ); 
//       for ( int l=0; l<indices.size(); l++)
//         std::cout << "K[ " << i << " , " << indices[l] << " ]  =  " << values[l] << "    ";
//       std::cout << endl;
//   }
//   
//   std::cout << "BBB" << std::endl;
//   for ( long int i=0; i<B->getRowMap()->getGlobalNumElements(); i++ )
//   {
//       B->getLocalRowView( i, indices, values ); 
//       for ( int l=0; l<indices.size(); l++)
//         std::cout << "K[ " << i << " , " << indices[l] << " ]  =  " << values[l] << "    ";
//       std::cout << endl;
//   }
  
  

    if ( firstTime_ )
    {
        komplex_->finalizeMatrix();
        firstTime_ = false;
    }
    
//     std::cout << "komplex" << std::endl;
//     Epetra_CrsMatrix AA = *komplex_->getMatrix();
//     AA.Scale( 1.0/16.0 );
//     std::cout << AA << std::endl;

    return true;
}
// =============================================================================
bool
Ginla::FDM::LocaSystem::Default::
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
Ginla::FDM::LocaSystem::Default::
getJacobian() const
{
    return komplex_->getMatrix();
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
Ginla::FDM::LocaSystem::Default::getPreconditioner() const
{
    return preconditioner_;
}
// =============================================================================
bool
Ginla::FDM::LocaSystem::Default::
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
Ginla::FDM::LocaSystem::Default::
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
Ginla::FDM::LocaSystem::Default::
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
Ginla::FDM::LocaSystem::Default::
releaseLocaStepper()
{
    stepper_ = Teuchos::null;
}
// =============================================================================
// function used by LOCA
void
Ginla::FDM::LocaSystem::Default::
printSolution ( const  Epetra_Vector &x,
                double conParam
              )
{
    // define vector
    const Teuchos::RCP<const Ginla::FDM::State> state =
        this->createFdmState_( x );

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
Ginla::FDM::LocaSystem::Default::
printSolution ( const Epetra_Vector & x,
                const std::string   & filenameAppendix
              ) const
{
    const Teuchos::RCP<const Ginla::FDM::State> state =
        this->createFdmState_( x );

    stateWriter_->write( state,
                         stepper_->getStepNumber(),
                         "-" + filenameAppendix
                       );
    return;
}
// =============================================================================
void
Ginla::FDM::LocaSystem::Default::
printSolutionOneParameterContinuation ( const Teuchos::RCP<const Ginla::FDM::State> & state )
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
Ginla::FDM::LocaSystem::Default::
printSolutionTurningPointContinuation ( const Teuchos::RCP<const Ginla::FDM::State> & state )
{
    static bool isSolution = false;

    // alternate between solution and nullvector
    isSolution = !isSolution;
  
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
Ginla::FDM::LocaSystem::Default::
writeContinuationStats ( const Teuchos::RCP<const Ginla::FDM::State> & state )
{   
    TEUCHOS_ASSERT( !statsWriter_.is_null() );
    
    Teuchos::RCP<Teuchos::ParameterList> paramList = statsWriter_->getListNonConst();

    TEUCHOS_ASSERT( !paramList.is_null() );
                    
    paramList->set( "0step", stepper_->getStepNumber() );
    paramList->set( "2#nonlinear steps", stepper_->getSolver()->getNumIterations() );
  
    // put the parameter list into statsWriter_
    std::string labelPrepend = "1";
    Ginla::Helpers::appendToTeuchosParameterList( *paramList,
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
Ginla::FDM::LocaSystem::Default::
setOutputDir ( const string & directory )
{
    stateWriter_->setOutputDir( directory );
    return;
}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
Ginla::FDM::LocaSystem::Default::
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
Ginla::FDM::LocaSystem::Default::
getMap() const
{
    return komplex_->getRealMap();
}
// =============================================================================
Teuchos::RCP<Ginla::State::Virtual>
Ginla::FDM::LocaSystem::Default::
createState( const Epetra_Vector & x
           ) const
{
    return this->createFdmState_( x );
}
// =============================================================================
Teuchos::RCP<Ginla::FDM::State>
Ginla::FDM::LocaSystem::Default::
createFdmState_( const Epetra_Vector & x
               ) const
{
    Teuchos::RCP<ComplexVector> psi = komplex_->real2complex ( x );
    return Teuchos::rcp( new Ginla::FDM::State( psi, glOperator_->getGrid() ) );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::FDM::LocaSystem::Default::
createSystemVector( const Ginla::State::Virtual & state
                  ) const
{
    return komplex_->complex2real( state.getPsi() );
}
// =============================================================================