/*
 * GL::LocaSystem::Bordered.cpp
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl\"omer
 */

#include "Ginla_LocaSystem_Bordered.h"

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

// =============================================================================
// Default constructor
Ginla::LocaSystem::Bordered::
Bordered ( const Teuchos::RCP<Ginla::Operator::Virtual>           & glOperator,
           const Teuchos::RCP<const Epetra_Comm>                  & eComm,
           const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > & complexMap,
           const Teuchos::RCP<Ginla::IO::StatsWriter>             & statsWriter,
           const Teuchos::RCP<Ginla::IO::StateWriter>             & stateWriter
         ) :
        glSystem_ ( glOperator,
                    eComm,
                    complexMap,
                    statsWriter,
                    stateWriter ),
        regularMap_ (  glSystem_.getMap() ),
        extendedMap_ ( createExtendedRealMap ( *regularMap_ ) ),
        jacobian_ ( new Epetra_CrsMatrix ( Copy, *extendedMap_, 0 ) ),
        firstTime_ ( true ),
        importFromExtendedMap_(  Epetra_Import(*regularMap_,*extendedMap_) ),
        importFromRegularMap_( Epetra_Import(*extendedMap_,*regularMap_) )
{
  return;
}
// =============================================================================
// Destructor
Ginla::LocaSystem::Bordered::
~Bordered()
{
}
// =============================================================================
bool
Ginla::LocaSystem::Bordered::
computeF ( const Epetra_Vector & x,
           Epetra_Vector       & FVec,
           const NOX::Epetra::Interface::Required::FillType fillFlag
         )
{
    TEST_FOR_EXCEPTION ( !regularMap_.is_valid_ptr() || regularMap_.is_null(),
                         std::logic_error,
                         "regularMap_ not properly initialized." );

    TEST_FOR_EXCEPTION ( !extendedMap_.is_valid_ptr() || extendedMap_.is_null(),
                         std::logic_error,
                         "extendedMap_ not properly initialized." );

    // make sure that the input and output vectors are correctly mapped
    TEST_FOR_EXCEPTION ( !x.Map().SameAs ( *extendedMap_ ),
                         std::logic_error,
                         "Maps of x and the computed real-valued map do not coincide. "
                         << "Check, for example, the number of elements "
                         << "(" << x.Map().NumGlobalElements() << " for x vs. "
                         << extendedMap_->NumGlobalElements() << " for extendedMap_)." );

    TEST_FOR_EXCEPTION ( !FVec.Map().SameAs ( *extendedMap_ ),
                         std::logic_error,
                         "Maps of FVec and the computed real-valued map do not coincide."
                         << "Check, for example, the number of elements "
                         << "(" << FVec.Map().NumGlobalElements() << " for FVec vs. "
                         << extendedMap_->NumGlobalElements() << " for extendedMap_)." );
    
    // strip off the phase constraint variable
    Epetra_Vector tmp ( *regularMap_ );
    tmp.Import( x, importFromExtendedMap_, Insert );

    double eta = x[ x.GlobalLength()-1 ];

    // compute GL
    Epetra_Vector shortFVec ( *regularMap_ );
    glSystem_.computeF ( tmp, shortFVec, fillFlag );

    // add -i\eta\psi
    Teuchos::RCP<Ginla::State> state = glSystem_.createState( tmp );
    state->getValuesNonConst()->scale( -IM*eta );
    TEUCHOS_ASSERT_EQUALITY( 0, tmp.Update( 1.0, *glSystem_.createSystemVector( *state ), 1.0 ) );
    
    // copy over and add phase condition
    FVec.Import( shortFVec, importFromRegularMap_, Insert );

    // set last entry
    FVec.ReplaceGlobalValue ( shortFVec.GlobalLength(), 0, 0.0 );

    return true;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
Ginla::LocaSystem::Bordered::
createExtendedRealMap ( const Epetra_BlockMap & realMap ) const
{
    // fill up realMapGIDs
    int numMyElements = realMap.NumMyElements();
    Teuchos::Array<int> myElements ( numMyElements );
    realMap.MyGlobalElements ( myElements.getRawPtr() );

    // Unconditionally put the phase constraint on the first process.
    int myPID = realMap.Comm().MyPID();
    if ( myPID==0 )
    {
        int n = realMap.NumGlobalElements();
        // extend the GIDs by the phase constraint
        myElements.append ( n );
    }

    int numGlobalElements = realMap.NumGlobalElements() + 1;
    return Teuchos::rcp ( new Epetra_Map ( numGlobalElements,
                                           myElements.length(),
                                           myElements.getRawPtr(),
                                           realMap.IndexBase(),
                                           realMap.Comm() )
                        );
}
// =============================================================================
bool
Ginla::LocaSystem::Bordered::
computeJacobian ( const Epetra_Vector   & x,
                        Epetra_Operator & Jac
                )
{
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    tmp.Import( x, importFromExtendedMap_, Insert );

    // Strip down Jac, too?
    // --   Not really necessary as it's not being used anyway.

    // compute the underlying Jacobian
    glSystem_.computeJacobian ( tmp, Jac );
    
    // compute the values of the Jacobian
    createJacobian ( x );

    return true;
}
// =============================================================================
bool
Ginla::LocaSystem::Bordered::
computePreconditioner ( const Epetra_Vector    & x,
                        Epetra_Operator        & Prec,
                        Teuchos::ParameterList * precParams )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Use explicit Jacobian only for this test problem!" );
    return true;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
Ginla::LocaSystem::Bordered::getJacobian() const
{
    return jacobian_;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
Ginla::LocaSystem::Bordered::getPreconditioner() const
{
    return preconditioner_;
}
// =============================================================================
// It also incorporates a phase condition.
void
Ginla::LocaSystem::Bordered::
createJacobian ( const Epetra_Vector & x )
{
//     const Ginla::Komplex & komplex = *glSystem_.getKomplex();
  
    TEST_FOR_EXCEPTION ( !extendedMap_.is_valid_ptr() || extendedMap_.is_null(),
                         std::logic_error,
                         "extendedMap_ not properly initialized." );

    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    tmp.Import( x, importFromExtendedMap_, Insert );
    
    double eta = x[ x.GlobalLength()-1 ];

    // TODO don't explicitly construct psi? get1dCopy on the rhs
    Teuchos::RCP<Ginla::State>              state   = glSystem_.createState( tmp );
//     Teuchos::ArrayRCP<const double_complex> psiView = psi->get1dView();

    // get the unbordered Jacobian
    Teuchos::RCP<Epetra_CrsMatrix> regularJacobian = glSystem_.getJacobian();
    
    // add -i*eta on the diagonal
    Epetra_Vector newDiag ( tmp );
    regularJacobian->ExtractDiagonalCopy( newDiag );
    Teuchos::RCP<Ginla::State> diag = glSystem_.createState( newDiag );
    diag->getValuesNonConst()->putScalar( -IM*eta );
    newDiag.Update( 1.0, *glSystem_.createSystemVector(*diag), 1.0 );
    regularJacobian->ReplaceDiagonalValues( newDiag );
    
    // TODO: Conversion to real-valued vector in one go?
    // right bordering: (phi:=) -i*psi
    Ginla::State phi = *state;
    phi.getValuesNonConst()->scale ( -IM );
    Teuchos::RCP<Epetra_Vector> rightBorder = glSystem_.createSystemVector( phi );

    // Get the lower bordering
    //     
    //    Im( psi_{old}^H, dpsi ) 
    //  = \Im( \sum_{k=1}^n \psi^*_k dpsi_k )
    //  = \sum_{k=1}^n ( -Im(\psi_k) Re(dpsi_k) + Re(\psi_k) Im(dpsi_k)
    //
    // Hence, the bordering tranlates into
    //
    //  lb = IM*psi_{old}
    //
    // in the complex state notation.
    Ginla::State lb_state = *state;
    lb_state.getValuesNonConst()->scale( IM );
    Teuchos::RCP<Epetra_Vector> lowerBorder = glSystem_.createSystemVector( lb_state );

    // corner element
    double d = 0.0;

    // create the bordered Jacobian out of this
    fillBorderedMatrix ( jacobian_,
                         regularJacobian,
                         *rightBorder,
                         *lowerBorder,
                         d,
                         firstTime_ );

    if ( firstTime_ )
    {
        TEUCHOS_ASSERT_EQUALITY ( 0, jacobian_->FillComplete() );
        TEUCHOS_ASSERT_EQUALITY ( 0, jacobian_->OptimizeStorage() );
        firstTime_ = false;
    }

    return;
}
// =============================================================================
bool
Ginla::LocaSystem::Bordered::
computeShiftedMatrix ( double alpha,
                       double beta,
                       const Epetra_Vector & x,
                       Epetra_Operator     & A )
{
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    tmp.Import( x, importFromExtendedMap_, Insert );
    
    // compute the underlying Jacobian
    glSystem_.computeJacobian ( tmp, A );
  
    // compute the values of the Jacobian
    createJacobian ( x );

    jacobian_->Scale ( alpha );
    //  jacobian_->FillComplete();

    Epetra_Vector newDiag ( x );
    Epetra_Vector unitVector ( x );
    unitVector.PutScalar ( 1.0 );
    //  newDiag.PutScalar(0.0);
    jacobian_->ExtractDiagonalCopy ( newDiag );
    newDiag.Update ( beta, unitVector, 1.0 );
    jacobian_->ReplaceDiagonalValues ( newDiag );

    return true;
}
// =============================================================================
// function used by LOCA
void
Ginla::LocaSystem::Bordered::setParameters ( const LOCA::ParameterVector &p )
{
    glSystem_.setParameters ( p );
}
// =============================================================================
void
Ginla::LocaSystem::Bordered::
setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper )
{
    glSystem_.setLocaStepper ( stepper );
}
// =============================================================================
void
Ginla::LocaSystem::Bordered::
releaseLocaStepper()
{
    glSystem_.releaseLocaStepper();
}
// =============================================================================
// function used by LOCA
void
Ginla::LocaSystem::Bordered::
printSolution ( const  Epetra_Vector &x,
                double conParam )
{
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    tmp.Import( x, importFromExtendedMap_, Insert );

    glSystem_.printSolution ( tmp, conParam );
}
// =============================================================================
void
Ginla::LocaSystem::Bordered::
printSolution ( const Epetra_Vector & x,
                const std::string   & filenameAppendix
              ) const
{
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    tmp.Import( x, importFromExtendedMap_, Insert );
    
    glSystem_.printSolution ( tmp, filenameAppendix );
    return;
}
// =============================================================================
// function used by LOCA
void Ginla::LocaSystem::Bordered::
setOutputDir ( const string &directory )
{
    glSystem_.setOutputDir ( directory );
}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
Ginla::LocaSystem::Bordered::
create_CommInt ( const Teuchos::RCP<const Epetra_Comm> &epetraComm )
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
Ginla::LocaSystem::Bordered::getMap() const
{
    return extendedMap_;
}
// =============================================================================
void
Ginla::LocaSystem::Bordered::
fillBorderedMatrix ( const Teuchos::RCP<      Epetra_CrsMatrix> & extendedMatrix,
                     const Teuchos::RCP<const Epetra_CrsMatrix> & regularMatrix,
                     const Epetra_Vector                        & rightBorder,
                     // TODO Declare the following const as soon as Trilinos allows (ReplaceGlobalValues)
                     Epetra_Vector                              & lowerBorder,
                     double                                       d,
                     bool                                         firstTime
                   ) const
{
    TEUCHOS_ASSERT ( regularMatrix.is_valid_ptr()  && !regularMatrix.is_null() );
    TEUCHOS_ASSERT ( extendedMatrix.is_valid_ptr() && !extendedMatrix.is_null() );

    int m = regularMatrix->NumGlobalRows();
    int n = regularMatrix->NumGlobalCols();

    // check if the sizes all match
    TEUCHOS_ASSERT_EQUALITY ( m+1, extendedMatrix->NumGlobalRows() );
    TEUCHOS_ASSERT_EQUALITY ( n+1, extendedMatrix->NumGlobalCols() );

    // make sure the maps coincide
    TEUCHOS_ASSERT ( lowerBorder.Map().SameAs ( regularMatrix->OperatorDomainMap() ) );
    TEUCHOS_ASSERT ( rightBorder.Map().SameAs ( regularMatrix->OperatorRangeMap() ) );

    int numMyRows = regularMatrix->NumMyRows();

    // fill the matrix with the entries
    int numRowNonZeros;

    int maxNumEntries = regularMatrix->MaxNumEntries() + 1; // count the last column in
    int    * indices = new int   [maxNumEntries];
    double * values  = new double[maxNumEntries];
    for ( int myRow=0; myRow<numMyRows; myRow++ )
    {
        // extract row view
        TEUCHOS_ASSERT_EQUALITY ( 0, regularMatrix->ExtractMyRowView ( myRow, numRowNonZeros, values, indices ) );

        // *indices are given in global indexing.
        int globalRow = extendedMatrix->Map().GID ( myRow );

        // Write the data to the new matrix.
        // Only panic for negative return codes.
        TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( extendedMatrix, globalRow, numRowNonZeros, values, indices, firstTime ) );

        // add last column
        double val = rightBorder[globalRow];
        TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( extendedMatrix, globalRow, 1, &val, &n, firstTime ) );
    }


    // set the last row
    // put the last row piece by piece
    int    numMyElements = lowerBorder.Map().NumMyElements();
    std::vector<double> myLowerBorderValues ( numMyElements );
    // TODO only use ExtractView
    lowerBorder.ExtractCopy ( & ( myLowerBorderValues[0] ) );
    int * myLowerBorderIndices = new int[numMyElements];
    for ( int k=0; k<numMyElements; k++ )
        myLowerBorderIndices[k] = lowerBorder.Map().GID ( k );
    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( extendedMatrix, n, numMyElements, & ( myLowerBorderValues[0] ), myLowerBorderIndices, firstTime ) );

    // set the last element d
    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( extendedMatrix, n, 1, &d, &n, firstTime ) );

    return;
}
// =============================================================================
int
Ginla::LocaSystem::Bordered::
PutRow ( const Teuchos::RCP<Epetra_CrsMatrix> A,
         const int                            globalRow,
         const int                            numIndices,
         double                             * values,
         int                                * indices,
         const bool                           firstTime
       ) const
{
    if ( firstTime )
        return A->InsertGlobalValues ( globalRow, numIndices, values, indices );
    else
        return A->ReplaceGlobalValues ( globalRow, numIndices, values, indices );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::LocaSystem::Bordered::
createSystemVector( const Ginla::State & state ) const
{ 
    Teuchos::RCP<Epetra_Vector> x =
        Teuchos::rcp( new Epetra_Vector( *extendedMap_ ) );

    x->Import( *glSystem_.createSystemVector( state ),
               importFromRegularMap_,
               Insert );

    // set last entry
    x->ReplaceGlobalValue ( x->GlobalLength()-1, 0, state.getChi() );
                
    return x;
}
// =============================================================================
Teuchos::RCP<Ginla::State>
Ginla::LocaSystem::Bordered::
createState( const Epetra_Vector & x ) const
{
  // TODO Use chi as well!
  
  Epetra_Vector tmp ( *regularMap_ );
  tmp.Import( x, importFromExtendedMap_, Insert );
  
  return glSystem_.createState( tmp );
}
// =============================================================================