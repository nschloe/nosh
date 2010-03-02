/*
 * GL::LinearSystem::Bordered.cpp
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl\"omer
 */

#include "GL_LinearSystem_Bordered.h"

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
GL::LinearSystem::Bordered::Bordered ( GinzburgLandau::GinzburgLandau &gl,
                                       const Teuchos::RCP<const Epetra_Comm> eComm,
                                       const Teuchos::RCP<const ComplexVector> psi,
                                       const std::string outputDir,
                                       const std::string outputDataFileName,
                                       const std::string outputFileFormat,
                                       const std::string solutionFileNameBase,
                                       const std::string nullvectorFileNameBase,
                                       const unsigned int maxStepNumberDecimals
                                     ) :
        glSystem_ ( gl, eComm, psi, outputDir, outputDataFileName, outputFileFormat,
                    solutionFileNameBase, nullvectorFileNameBase, maxStepNumberDecimals ),
        regularMap_ (  glSystem_.getRealMap() ),
        extendedMap_ ( createExtendedRealMap ( *regularMap_ ) ),
        jacobian_ ( new Epetra_CrsMatrix ( Copy, *extendedMap_, 0 ) ),
        solution_ ( new Epetra_Vector ( *extendedMap_ ) ),
        maxStepNumberDecimals_ ( maxStepNumberDecimals ),
        firstTime_ ( true )
{
    Teuchos::RCP<Epetra_Vector> tmp = glSystem_.getGlKomplex()->complex2real ( *psi );

    for ( int k=0; k<tmp->MyLength(); k++ )
        (*solution_)[k] = ( *tmp ) [tmp->Map().GID ( k ) ];

    int n = solution_->GlobalLength();
    solution_->ReplaceGlobalValue ( n-1, 0, 0.0 );

    // Initialize the format for the the continuation step number.
    // Here: 00012 for step no. 12, if maxStepNumberDecimals_=5.
//     stepNumFileNameFormat_ = boost::str ( boost::format ( "%%|0%d|" ) % maxStepNumberDecimals_ );
}
// =============================================================================
// Destructor
GL::LinearSystem::Bordered::~Bordered()
{
}
// =============================================================================
bool
GL::LinearSystem::Bordered::computeF ( const Epetra_Vector & x,
                                       Epetra_Vector       & FVec,
                                       const NOX::Epetra::Interface::Required::FillType fillFlag )
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
    
    // TODO replace by {im,ex}porter
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    for ( int k=0; k<tmp.MyLength(); k++ )
        tmp[k] = x[x.Map().GID ( k ) ];

    Epetra_Vector shortFVec ( *regularMap_ );
    glSystem_.computeF ( tmp, shortFVec, fillFlag );

    // copy over and add phase condition
    // TODO replace by {im,ex}porter
    for ( int k=0; k<shortFVec.MyLength(); k++ )
        FVec[k] = shortFVec[shortFVec.Map().GID ( k ) ];

    FVec.ReplaceGlobalValue ( shortFVec.GlobalLength(), 0, 0.0 );

    return true;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
GL::LinearSystem::Bordered::createExtendedRealMap ( const Epetra_BlockMap & realMap ) const
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
GL::LinearSystem::Bordered::computeJacobian ( const Epetra_Vector & x,
                                          Epetra_Operator     & Jac
                                        )
{
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    for ( int k=0; k<tmp.MyLength(); k++ )
        tmp[k] = x[x.Map().GID ( k ) ];

    // TODO Strip down Jac, too.
    // --   Not really necessary as it's not being used anyway.

    // compute the underlying Jacobian
    glSystem_.computeJacobian ( tmp, Jac );
    
    // compute the values of the Jacobian
    createJacobian ( x );

    return true;
}
// =============================================================================
bool
GL::LinearSystem::Bordered::computePreconditioner ( const Epetra_Vector    & x,
                                                Epetra_Operator        & Prec,
                                                Teuchos::ParameterList * precParams )
{
    TEST_FOR_EXCEPTION ( true,
                         std::logic_error,
                         "Use explicit Jacobian only for this test problem!" );
    return true;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
GL::LinearSystem::Bordered::getSolution() const
{
    return solution_;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GL::LinearSystem::Bordered::getJacobian() const
{
    return jacobian_;
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
GL::LinearSystem::Bordered::getPreconditioner() const
{
    return preconditioner_;
}
// =============================================================================
// It also incorporates a phase condition.
void
GL::LinearSystem::Bordered::createJacobian ( const Epetra_Vector & x )
{
    TEST_FOR_EXCEPTION ( !extendedMap_.is_valid_ptr() || extendedMap_.is_null(),
                         std::logic_error,
                         "extendedMap_ not properly initialized." );

    // TODO replace by {im,ex}porter
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    for ( int k=0; k<tmp.MyLength(); k++ )
        tmp[k] = x[x.Map().GID ( k ) ];

    // TODO don't explicitly construct psi? get1dCopy on the rhs
    Teuchos::RCP<ComplexVector>             psi     = glSystem_.getGlKomplex()->real2complex ( tmp );
    Teuchos::ArrayRCP<const double_complex> psiView = psi->get1dView();

    // get the unbordered Jacobian
    Teuchos::RCP<const Epetra_CrsMatrix> regularJacobian = glSystem_.getJacobian();
    
    // TODO: Conversion to real-valued vector in one go?
    // right bordering: (phi:=) -i*psi
    ComplexVector phi = *psi;
    phi.scale ( std::complex<double> ( 0.0,-1.0 ) );
    Teuchos::RCP<Epetra_Vector> rightBorder =
        glSystem_.getGlKomplex()->complex2real ( phi );

    // Get the lower bordering  Im( psi_{old}^H, dpsi ).
    Teuchos::RCP<Epetra_Vector> lowerBorder =
        glSystem_.getGlKomplex()->imagScalarProductCoeff ( psi );

    // corner element
    double d = 0.0;

    // create the bordered Jacobian out of this
    fillBorderedMatrix ( jacobian_, regularJacobian, *rightBorder, *lowerBorder, d, firstTime_ );

    if ( firstTime_ )
    {
        TEST_FOR_EXCEPT ( 0 != jacobian_->FillComplete() );
        TEST_FOR_EXCEPT ( 0 != jacobian_->OptimizeStorage() );
        firstTime_ = false;
    }

    return;
}
// =============================================================================
bool
GL::LinearSystem::Bordered::computeShiftedMatrix ( double alpha,
                                               double beta,
                                               const Epetra_Vector   & x,
                                                     Epetra_Operator & A )
{
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    for ( int k=0; k<tmp.MyLength(); k++ )
        tmp[k] = x[x.Map().GID ( k ) ];
    
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
GL::LinearSystem::Bordered::setParameters ( const LOCA::ParameterVector &p )
{
    glSystem_.setParameters ( p );
}
// =============================================================================
void
GL::LinearSystem::Bordered::setLocaStepper ( const Teuchos::RCP<const LOCA::Stepper> stepper )
{
    glSystem_.setLocaStepper ( stepper );
}
// =============================================================================
void
GL::LinearSystem::Bordered::releaseLocaStepper()
{
    glSystem_.releaseLocaStepper();
}
// =============================================================================
// function used by LOCA
void
GL::LinearSystem::Bordered::printSolution ( const  Epetra_Vector &x,
                                        double conParam )
{
    // TODO replace by {im,ex}porter
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    for ( int k=0; k<tmp.MyLength(); k++ )
        tmp[k] = x[x.Map().GID ( k ) ];

    glSystem_.printSolution ( tmp, conParam );
}
// =============================================================================
// function used by LOCA
void GL::LinearSystem::Bordered::setOutputDir ( const string &directory )
{
    glSystem_.setOutputDir ( directory );
}
// =============================================================================
void
GL::LinearSystem::Bordered::writeSolutionToFile ( const Epetra_Vector & x,
                                              const std::string   & filePath
                                            ) const
{
    // TODO replace by {im,ex}porter
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    for ( int k=0; k<tmp.MyLength(); k++ )
        tmp[k] = x[x.Map().GID ( k ) ];

    glSystem_.writeSolutionToFile ( tmp, filePath );
}
// =============================================================================
void
GL::LinearSystem::Bordered::writeAbstractStateToFile ( const Epetra_Vector & x,
                                                   const std::string   & filePath
                                                 ) const
{
    // TODO replace by {im,ex}porter
    // strip off the phase constraint
    Epetra_Vector tmp ( *regularMap_ );
    for ( int k=0; k<tmp.MyLength(); k++ )
        tmp[k] = x[x.Map().GID ( k ) ];

    glSystem_.writeAbstractStateToFile ( tmp, filePath );
}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
GL::LinearSystem::Bordered::create_CommInt ( const Teuchos::RCP<const Epetra_Comm> &epetraComm )
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
// TODO delete?
const Teuchos::RCP<const GL::Komplex>
GL::LinearSystem::Bordered::getGlKomplex() const
{
    return glSystem_.getGlKomplex();
}
// =============================================================================
double
GL::LinearSystem::Bordered::getH0() const
{
    return glSystem_.getH0();
}
// =============================================================================
const Teuchos::RCP<const Epetra_Map>
GL::LinearSystem::Bordered::getMap() const
{
    return extendedMap_;
}
// =============================================================================
void
GL::LinearSystem::Bordered::setH0 ( const double h0 )
{
    glSystem_.setH0 ( h0 );
}
// =============================================================================
void
GL::LinearSystem::Bordered::setScaling ( const double scaling )
{
    glSystem_.setScaling ( scaling );
}
// =============================================================================
void
GL::LinearSystem::Bordered::setChi ( const double chi )
{
    glSystem_.setChi ( chi );
}
// =============================================================================
void
GL::LinearSystem::Bordered::fillBorderedMatrix ( const Teuchos::RCP<      Epetra_CrsMatrix> & extendedMatrix,
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
GL::LinearSystem::Bordered::PutRow ( const Teuchos::RCP<Epetra_CrsMatrix> A,
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
