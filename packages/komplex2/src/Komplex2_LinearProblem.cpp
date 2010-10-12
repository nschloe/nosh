/*
 * GL::Komplex.cpp
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl\"omer
 */

#include "Komplex2_LinearProblem.h"
#include "Komplex2_DoubleMatrix.h"

#include <Epetra_Map.h>
#include <Teuchos_DefaultComm.hpp> // for Teuchos::SerialComm
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif // HAVE_MPI

// =============================================================================
// Default constructor
Komplex2::LinearProblem::
LinearProblem ( const Teuchos::RCP<const Epetra_Comm>                  eComm,
                const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > complexRowMap
              ) :
        EComm_ ( eComm ),
        TComm_ ( this->create_CommInt_ ( eComm ) ),
        RealMap_ ( this->createRealMap_ ( complexRowMap ) ),
        ComplexMap_ ( complexRowMap ),
        realMatrix_ ( Teuchos::rcp( new Epetra_CrsMatrix (Copy,*RealMap_,0) ) )
{
}
// =============================================================================
Komplex2::LinearProblem::
~LinearProblem()
{
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
Komplex2::LinearProblem::
getRealMap() const
{
    return RealMap_;
}
// =============================================================================
Teuchos::RCP<const ComplexMap>
Komplex2::LinearProblem::
getComplexMap() const
{
    return ComplexMap_;
}
// =============================================================================
// converts a real-valued vector to a complex-valued psi vector
Teuchos::RCP<ComplexVector>
Komplex2::LinearProblem::
real2complex ( const Epetra_Vector & x ) const
{
    TEST_FOR_EXCEPTION ( !RealMap_.is_valid_ptr() || RealMap_.is_null(),
                         std::logic_error,
                         "RealMap_ not properly initialized." );

    TEST_FOR_EXCEPTION ( !x.Map().SameAs ( *RealMap_ ),
                         std::logic_error,
                         "Maps for real-valued vectors do not coincide. "
                         << "Check, for example, the number of elements "
                         << "(" << x.Map().NumGlobalElements() << " for x vs. "
                         << RealMap_->NumGlobalElements() << " for RealMap_)." );

    TEST_FOR_EXCEPTION ( !ComplexMap_.is_valid_ptr() || ComplexMap_.is_null(),
                         std::logic_error,
                         "ComplexMap_ has not been properly initialized." );

    Teuchos::RCP<ComplexVector> z = Teuchos::rcp ( new ComplexVector ( ComplexMap_ ) );
    Teuchos::ArrayRCP<double_complex> zView = z->get1dViewNonConst();
    for ( unsigned int k=0; k < z->getLocalLength(); k++ )
        zView[k] = double_complex ( x[2*k], x[2*k+1] );

    return z;
}
// =============================================================================
// converts a real-valued vector to a complex-valued psi vector
// Teuchos::RCP<ComplexMultiVector>
// GL::Komplex::real2complex ( const Epetra_MultiVector & x ) const
// {
//     TEST_FOR_EXCEPTION ( !RealMap_.is_valid_ptr() || RealMap_.is_null(),
//                          std::logic_error,
//                          "RealMap_ not properly initialized." );
// 
//     TEST_FOR_EXCEPTION ( !x.Map().SameAs ( *RealMap_ ),
//                          std::logic_error,
//                          "Maps for real-valued vectors do not coincide. "
//                          << "Check, for example, the number of elements "
//                          << "(" << x.Map().NumGlobalElements() << " for x vs. "
//                          << RealMap_->NumGlobalElements() << " for RealMap_)." );
// 
//     TEST_FOR_EXCEPTION ( !ComplexMap_.is_valid_ptr() || ComplexMap_.is_null(),
//                          std::logic_error,
//                          "ComplexMap_ has not been properly initialized." );
// 
//     Teuchos::RCP<ComplexVector> z = Teuchos::rcp ( new ComplexVector ( ComplexMap_ ) );
//     Teuchos::ArrayRCP<double_complex> zView = z->get1dViewNonConst();
//     for ( unsigned int k=0; k < z->getLocalLength(); k++ )
//         zView[k] = double_complex ( x[2*k], x[2*k+1] );
// 
//     return z;
// }
// =============================================================================
// converts a real-valued vector to a complex-valued psi vector
Teuchos::RCP<Epetra_Vector>
Komplex2::LinearProblem::
complex2real ( const ComplexVector & complexVec ) const
{
    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp ( new Epetra_Vector ( *RealMap_ ) );
    this->complex2real( complexVec, *x );
    return x;
}
// =============================================================================
// converts a real-valued vector to a complex-valued psi vector
void
Komplex2::LinearProblem::
complex2real ( const ComplexVector & complexVec,
                     Epetra_Vector & x
             ) const
{
    TEUCHOS_ASSERT ( ComplexMap_.is_valid_ptr()
                     && !ComplexMap_.is_null()
                     && complexVec.getMap()->isSameAs ( *ComplexMap_ ) );

    TEUCHOS_ASSERT( !RealMap_.is_null()
                    && x.Map().SameAs( *RealMap_ ) );

    Teuchos::ArrayRCP<const double_complex> complexVecView =
        complexVec.get1dView();
    for ( unsigned int k = 0; k < ComplexMap_->getNodeNumElements(); k++ )
    {
        x[2*k]   = std::real ( complexVecView[k] );
        x[2*k+1] = std::imag ( complexVecView[k] );
    }
    return;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Komplex2::LinearProblem::
complex2real ( const Teuchos::RCP<const ComplexVector> & complexVecPtr ) const
{
  
    TEUCHOS_ASSERT ( complexVecPtr.is_valid_ptr() && !complexVecPtr.is_null() );
    return this->complex2real ( *complexVecPtr );
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
Komplex2::LinearProblem::
createRealMap_ ( const Teuchos::RCP<const ComplexMap> & ComplexMap
               ) const
{
    TEST_FOR_EXCEPTION ( !ComplexMap.is_valid_ptr() || ComplexMap.is_null(),
                         std::logic_error,
                         "ComplexMap has not been properly initialized." );

    // get view for the global indices of the global elements
    Teuchos::ArrayView<const Thyra::Ordinal> myComplexGIDs =
        ComplexMap->getNodeElementList();

    // Construct the map in such a way that all complex entries on processor K
    // are split up into real and imaginary part, which will both reside on
    // processor K again.
    unsigned int numMyComplexElements = myComplexGIDs.size();
    unsigned int numMyRealElements    = 2*numMyComplexElements;
    Teuchos::ArrayRCP<int> myRealGIDs ( numMyRealElements );
    for ( unsigned int k = 0; k < numMyComplexElements; k++ )
    {
        myRealGIDs[2*k  ] = 2 * myComplexGIDs[k];
        myRealGIDs[2*k+1] = 2 * myComplexGIDs[k] + 1;
    }
    
    return Teuchos::rcp ( new Epetra_Map ( -1,
                                           myRealGIDs.size(),
                                           myRealGIDs.getRawPtr(),
                                           ComplexMap->getIndexBase(),
                                           *EComm_
                                         )
                        );
}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
Komplex2::LinearProblem::
create_CommInt_ ( const Teuchos::RCP<const Epetra_Comm> &epetraComm )
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
        RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm =
            Teuchos::opaqueWrapper ( mpiEpetraComm->Comm() );
        set_extra_data ( mpiEpetraComm, "mpiEpetraComm", Teuchos::inOutArg ( rawMpiComm ) );
        RCP<const Teuchos::MpiComm<int> > mpiComm =
            rcp ( new Teuchos::MpiComm<int> ( rawMpiComm ) );
        return mpiComm;
    }
#else
    Teuchos::RCP<const Epetra_SerialComm>
    serialEpetraComm = Teuchos::rcp_dynamic_cast<const Epetra_SerialComm> ( epetraComm );
    if ( serialEpetraComm.get() )
    {
        Teuchos::RCP<const Teuchos::SerialComm<int> > serialComm =
            Teuchos::rcp ( new Teuchos::SerialComm<int>() );
        set_extra_data ( serialEpetraComm, "serialEpetraComm", Teuchos::inOutArg ( serialComm ) );
        return serialComm;
    }
#endif // HAVE_MPI

    // If you get here then the conversion failed!
    return Teuchos::null;
}
// =============================================================================
void
Komplex2::LinearProblem::
finalizeMatrix()
{
    TEUCHOS_ASSERT_EQUALITY ( 0, realMatrix_->FillComplete() );
    return;
}
// =============================================================================
void
Komplex2::LinearProblem::
zeroOutMatrix()
{
    TEUCHOS_ASSERT_EQUALITY ( 0, realMatrix_->PutScalar ( 0.0 ) );
}
// =============================================================================
void
Komplex2::LinearProblem::
update ( const Teuchos::RCP<const Komplex2::DoubleMatrix> AB,
         bool firstTime )
{
//     TEST_FOR_EXCEPTION( true,
//                         std::logic_error,
//                         "Function needs to be checked for mistakes -- there may be bugs."
//                       );

    int numMyElements = ComplexMap_->getNodeNumElements();
    Teuchos::ArrayRCP<const Thyra::Ordinal> indicesA;
    Teuchos::ArrayRCP<const double_complex> valuesA;
    Teuchos::ArrayRCP<const Thyra::Ordinal> indicesB;
    Teuchos::ArrayRCP<const double_complex> valuesB;

    for ( int row = 0; row < numMyElements; row++ )
    {
         Thyra::Ordinal globalRow = ComplexMap_->getGlobalElement ( row );

         AB->getMatrixA()->getLocalRowView( row,
                                             indicesA,
                                             valuesA
                                           );

         AB->getMatrixB()->getLocalRowView( row,
                                             indicesB,
                                             valuesB
                                           );

        this->updateGlobalRow ( globalRow,
                                indicesA(), valuesA(),
                                indicesB(), valuesB(),
                                firstTime
                              ); 
    }

    return;
}
// =============================================================================
void
Komplex2::LinearProblem::
updateGlobalRowA ( const unsigned int                                globalRow,
                   const Teuchos::ArrayView<const Thyra::Ordinal>  & indicesA,
                   const Teuchos::ArrayView<const double_complex>  & valuesA,
                   const bool                                        firstTime
                 )
{
    TEUCHOS_ASSERT ( realMatrix_.is_valid_ptr() && !realMatrix_.is_null() );
    TEUCHOS_ASSERT_INEQUALITY ( globalRow, <, ComplexMap_->getGlobalNumElements() );

    int numEntries = indicesA.size();
    Teuchos::ArrayRCP<int>    indicesReal( numEntries );
    Teuchos::ArrayRCP<double> valuesReal( numEntries );

    int realRow = 2*globalRow;
    // -------------------------------------------------------------------
    // insert the coefficients Re(A) of Re(z)
    for ( int k = 0; k < numEntries; k++ )
    {
        indicesReal[k] = 2 * indicesA[k];
        valuesReal[k]  = std::real ( valuesA[k] );
    }

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutGlobalRow_ ( realRow,
                                                             valuesReal(),
                                                             indicesReal(),
                                                             firstTime
                                                           )
                              );
    // -------------------------------------------------------------------
    // insert the coefficients -Im(A) of Im(z)    
    for ( int k = 0; k < numEntries; k++ )
    {
        indicesReal[k] = 2 * indicesA[k] + 1;
        valuesReal[k]  = -std::imag ( valuesA[k] );
    }

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutGlobalRow_ ( realRow,
                                                             valuesReal(),
                                                             indicesReal(),
                                                             firstTime
                                                           )
                              );
    // -------------------------------------------------------------------

    int imagRow = 2*globalRow+1;
    // -------------------------------------------------------------------
    // insert the coefficients Im(A) of Re(z)
    for ( int k = 0; k < numEntries; k++ )
    {
        indicesReal[k] = 2 * indicesA[k];
        valuesReal[k]  = std::imag ( valuesA[k] );
    }

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutGlobalRow_ ( imagRow,
                                                             valuesReal(),
                                                             indicesReal(),
                                                             firstTime
                                                           )
                              );
    // -------------------------------------------------------------------
    // insert the coefficients Re(A) of Im(z)    
    for ( int k = 0; k < numEntries; k++ )
    {
        indicesReal[k] = 2 * indicesA[k] + 1;
        valuesReal[k]  = std::real ( valuesA[k] );
    }

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutGlobalRow_ ( imagRow,
                                                             valuesReal(),
                                                             indicesReal(),
                                                             firstTime
                                                           )
                              );
    // -------------------------------------------------------------------
    
    return;
}
// =============================================================================
void
Komplex2::LinearProblem::
updateGlobalRowB ( const unsigned int                                globalRow,
                   const Teuchos::ArrayView<const Thyra::Ordinal>  & indicesB,
                   const Teuchos::ArrayView<const double_complex>  & valuesB,
                   const bool                                        firstTime
                 )
{
    TEUCHOS_ASSERT ( realMatrix_.is_valid_ptr() && !realMatrix_.is_null() );
    TEUCHOS_ASSERT_INEQUALITY ( globalRow, <, ComplexMap_->getGlobalNumElements() );

    int numEntries = indicesB.size();
    Teuchos::ArrayRCP<int>    indicesReal( numEntries );
    Teuchos::ArrayRCP<double> valuesReal( numEntries );

    int realRow = 2*globalRow;
    // -------------------------------------------------------------------
    // insert the coefficients Re(B) of Re(z)
    for ( int k = 0; k < numEntries; k++ )
    {
        indicesReal[k] = 2 * indicesB[k];
        valuesReal[k]  = std::real ( valuesB[k] );
    }

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutGlobalRow_ ( realRow,
                                                             valuesReal(),
                                                             indicesReal(),
                                                             firstTime
                                                           )
                              );
    // -------------------------------------------------------------------
    // insert the coefficients Im(B) of Im(z)
    for ( int k = 0; k < numEntries; k++ )
    {
        indicesReal[k] = 2 * indicesB[k] + 1;
        valuesReal[k]  = std::imag ( valuesB[k] );
    }

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutGlobalRow_ ( realRow, 
                                                             valuesReal(),
                                                             indicesReal(),
                                                             firstTime
                                                           )
                              );
    // -------------------------------------------------------------------

    int imagRow = 2*globalRow+1;
    // -------------------------------------------------------------------
    // insert the coefficients Im(B) of Re(z)
    for ( int k = 0; k < numEntries; k++ )
    {
        indicesReal[k] = 2 * indicesB[k];
        valuesReal[k]  = std::imag ( valuesB[k] );
    }

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutGlobalRow_ ( imagRow,
                                                             valuesReal(),
                                                             indicesReal(),
                                                             firstTime
                                                           )
                              );
    // -------------------------------------------------------------------
    // insert the coefficients -Re(B) of Im(z)
    for ( int k = 0; k < numEntries; k++ )
    {
        indicesReal[k] = 2 * indicesB[k] + 1;
        valuesReal[k]  = -std::real ( valuesB[k] );
    }

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, this->PutGlobalRow_ ( imagRow,
                                                             valuesReal(),
                                                             indicesReal(),
                                                             firstTime
                                                           )
                              );
    // -------------------------------------------------------------------

    return;
}
// =============================================================================
void
Komplex2::LinearProblem::
updateGlobalRow ( const unsigned int                                globalRow,
                  const Teuchos::ArrayView<const Thyra::Ordinal>  & indicesA,
                  const Teuchos::ArrayView<const double_complex>  & valuesA,
                  const Teuchos::ArrayView<const Thyra::Ordinal>  & indicesB,
                  const Teuchos::ArrayView<const double_complex>  & valuesB,
                  const bool                                        firstTime
                )
{
    this->updateGlobalRowA( globalRow,
                            indicesA, valuesA,
                            firstTime
                          );

    this->updateGlobalRowB( globalRow,
                            indicesB, valuesB,
                            firstTime
                          );

    return;
}
// =============================================================================
int
Komplex2::LinearProblem::
PutGlobalRow_ ( int                                    globalRow,
                const Teuchos::ArrayView<const double> values,
                const Teuchos::ArrayView<const int>    globalIndices,
                bool                                   firstTime
              )
{ 
    // We *need* to use the global variant here as the column map
    // of the real matrix is only determined at the first FillComplete(),
    // so InsertLocalValues() does not have any meaning.
    // This would usually work on one core, but inspecting the matrix
    // (,e.g., by Print()), shows invalid column indices.
    int numIndices = values.size();
    if ( firstTime )
        return realMatrix_->InsertGlobalValues ( globalRow,
                                                 numIndices,
                                                 values.getRawPtr(),
                                                 globalIndices.getRawPtr()
                                               );
    else
        return realMatrix_->SumIntoGlobalValues ( globalRow,
                                                  numIndices,
                                                  values.getRawPtr(),
                                                  globalIndices.getRawPtr()
                                                );
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
Komplex2::LinearProblem::
getMatrix() const
{
    return realMatrix_;
}
// =============================================================================
//Teuchos::RCP<Epetra_Vector>
//GL::Komplex::getRealCoefficients( const Teuchos::RCP<const ComplexVector> a,
//                                      const Teuchos::RCP<const ComplexVector> b
//                                    ) const
//{
//      TEUCHOS_ASSERT( a.is_valid_ptr() && !a.is_null() );
//      TEUCHOS_ASSERT( b.is_valid_ptr() && !b.is_null() );
//      TEUCHOS_ASSERT( a->getMap()->isSameAs(*ComplexMap_) );
//      TEUCHOS_ASSERT( b->getMap()->isSameAs(*ComplexMap_) );
//
//      Teuchos::RCP<Epetra_Vector> compound = Teuchos::rcp( new Epetra_Vector(*RealMap_) );
//
//      // TODO Handle without 1dViews
//      Teuchos::ArrayRCP<const double_complex> aView = a->get1dView();
//      Teuchos::ArrayRCP<const double_complex> bView = b->get1dView();
//      for ( unsigned int k=0; k<ComplexMap_->getNodeNumElements(); k++) {
//              int K = ComplexMap_->getGlobalElement(k);
//              compound->ReplaceMyValue( 2*k  , 0,  std::real(aView[K]) + std::real(bView[K]) );
//              compound->ReplaceMyValue( 2*k+1, 0, -std::imag(aView[K]) + std::imag(bView[K]) );
//      }
//
//      return compound;
//}
// =============================================================================
// Teuchos::RCP<Epetra_Vector>
// Komplex2::LinearProblem::
// imagScalarProductCoeff ( const Teuchos::RCP<const ComplexVector> a ) const
// {
//     TEUCHOS_ASSERT ( a.is_valid_ptr() && !a.is_null() );
//     TEUCHOS_ASSERT ( a->getMap()->isSameAs ( *ComplexMap_ ) );
// 
//     Teuchos::RCP<Epetra_Vector> compound = Teuchos::rcp ( new Epetra_Vector ( *RealMap_ ) );
// 
//     Teuchos::ArrayRCP<const double_complex> aView = a->get1dView();
//     for ( unsigned int k=0; k<ComplexMap_->getNodeNumElements(); k++ )
//     {
//         ( *compound ) [2*k]   = -std::imag ( aView[k] );
//         ( *compound ) [2*k+1] =  std::real ( aView[k] );
//     }
// 
//     return compound;
// }
// =============================================================================
