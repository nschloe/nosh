/*
 * GL::Komplex.cpp
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl\"omer
 */

#include "Ginla_Komplex.h"

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
Ginla::Komplex::
Komplex ( const Teuchos::RCP<const Epetra_Comm>                  eComm,
          const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > complexMap ) :
        EComm_ ( eComm ),
        TComm_ ( create_CommInt ( eComm ) ),
        RealMap_ ( createRealMap ( complexMap ) ),
        ComplexMap_ ( complexMap ),
        realMatrix_ ( new Epetra_CrsMatrix (Copy,*RealMap_,0) )
{
}
// =============================================================================
Ginla::Komplex::~Komplex()
{
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
Ginla::Komplex::getRealMap() const
{
    return RealMap_;
}
// =============================================================================
Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> >
Ginla::Komplex::getComplexMap() const
{
    return ComplexMap_;
}
// =============================================================================
// converts a real-valued vector to a complex-valued psi vector
Teuchos::RCP<ComplexVector>
Ginla::Komplex::real2complex ( const Epetra_Vector & x ) const
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
Ginla::Komplex::
complex2real ( const ComplexVector & complexVec ) const
{
    TEUCHOS_ASSERT ( ComplexMap_.is_valid_ptr()
                     && !ComplexMap_.is_null()
                     && complexVec.getMap()->isSameAs ( *ComplexMap_ ) );

    Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp ( new Epetra_Vector ( *RealMap_ ) );

    Teuchos::ArrayRCP<const double_complex> complexVecView = complexVec.get1dView();

    for ( unsigned int k = 0; k < ComplexMap_->getNodeNumElements(); k++ )
    {
        ( *x ) [2*k]   = std::real ( complexVecView[k] );
        ( *x ) [2*k+1] = std::imag ( complexVecView[k] );
    }
    return x;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::Komplex::
complex2real ( const Teuchos::RCP<const ComplexVector> & complexVecPtr ) const
{
  
    TEUCHOS_ASSERT ( complexVecPtr.is_valid_ptr() && !complexVecPtr.is_null() );
    return complex2real ( *complexVecPtr );
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
Ginla::Komplex::
createRealMap ( const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > & ComplexMap ) const
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
    Teuchos::Array<int> myRealGIDs ( numMyRealElements );
    for ( unsigned int k = 0; k < numMyComplexElements; k++ )
    {
        myRealGIDs[2*k  ] = 2 * myComplexGIDs[k];
        myRealGIDs[2*k+1] = 2 * myComplexGIDs[k] + 1;
    }

    return Teuchos::rcp ( new Epetra_Map ( numMyRealElements,
                                           myRealGIDs.size(),
                                           myRealGIDs.getRawPtr(),
                                           ComplexMap->getIndexBase(),
                                           *EComm_ )
                        );
}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
Ginla::Komplex::
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
Ginla::Komplex::finalizeMatrix()
{
    TEUCHOS_ASSERT_EQUALITY ( 0, realMatrix_->FillComplete() );
    TEUCHOS_ASSERT_EQUALITY ( 0, realMatrix_->OptimizeStorage() );
}
// =============================================================================
void
Ginla::Komplex::zeroOutMatrix()
{
    TEUCHOS_ASSERT_EQUALITY ( 0, realMatrix_->PutScalar ( 0.0 ) );
}
// =============================================================================
void
Ginla::Komplex::
updateRow ( const unsigned int                      row,
            const Teuchos::Array<int>             & indicesA,
            const Teuchos::Array<double_complex>  & valuesA,
            const Teuchos::Array<int>             & indicesB,
            const Teuchos::Array<double_complex>  & valuesB,
            const bool                              firstTime
          )
{
    TEUCHOS_ASSERT ( realMatrix_.is_valid_ptr() && !realMatrix_.is_null() );
    TEUCHOS_ASSERT_INEQUALITY ( row, <, ComplexMap_->getGlobalNumElements() );

    int numEntries;
    int * indicesAReal;
    double * valuesAReal;
    int * indicesAImag;
    double * valuesAImag;
    int * indicesBReal;
    double * valuesBReal;
    int * indicesBImag;
    double * valuesBImag;

    int realRow = 2*row;
    // -------------------------------------------------------------------
    // insert the coefficients Re(A) of Re(z)
    numEntries = indicesA.size();
    indicesAReal = new int[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        indicesAReal[k] = 2 * indicesA[k];

    valuesAReal = new double[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        valuesAReal[k] = std::real ( valuesA[k] );

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( realRow, numEntries, valuesAReal, indicesAReal, firstTime ) );
    delete[] valuesAReal;
    valuesAReal = NULL;
    delete[] indicesAReal;
    indicesAReal = NULL;
    // -------------------------------------------------------------------
    // insert the coefficients Re(B) of Re(z)
    numEntries = indicesB.size();
    indicesBReal = new int[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        indicesBReal[k] = 2 * indicesB[k];
    valuesBReal = new double[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        valuesBReal[k] = std::real ( valuesB[k] );

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( realRow, numEntries, valuesBReal, indicesBReal, firstTime ) );

    delete[] valuesBReal;
    valuesBReal = NULL;
    delete[] indicesBReal;
    indicesBReal = NULL;
    // -------------------------------------------------------------------
    // insert the coefficients -Im(A) of Im(z)
    numEntries = indicesA.size();
    indicesAImag = new int[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        indicesAImag[k] = 2 * indicesA[k] + 1;
    valuesAImag = new double[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        valuesAImag[k] = -std::imag ( valuesA[k] );

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( realRow, numEntries, valuesAImag, indicesAImag, firstTime ) );

    delete[] valuesAImag;
    valuesAImag = NULL;
    delete[] indicesAImag;
    indicesAImag = NULL;
    // -------------------------------------------------------------------
    // insert the coefficients Im(B) of Im(z)
    numEntries = indicesB.size();
    indicesBImag = new int[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        indicesBImag[k] = 2 * indicesB[k] + 1;
    valuesBImag = new double[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        valuesBImag[k] = std::imag ( valuesB[k] );

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( realRow, numEntries, valuesBImag, indicesBImag, firstTime ) );

    delete[] valuesBImag;
    valuesBImag = NULL;
    delete[] indicesBImag;
    indicesBImag = NULL;
    // -------------------------------------------------------------------


    int imagRow = 2*row+1;
    // -------------------------------------------------------------------
    // insert the coefficients Im(A) of Re(z)
    numEntries = indicesA.size();
    indicesAReal = new int[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        indicesAReal[k] = 2 * indicesA[k];
    valuesAImag = new double[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        valuesAImag[k] = std::imag ( valuesA[k] );

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( imagRow, numEntries, valuesAImag, indicesAReal, firstTime ) );

    delete[] valuesAImag;
    valuesAImag = NULL;
    delete[] indicesAReal;
    indicesAReal = NULL;
    // -------------------------------------------------------------------
    // insert the coefficients Im(B) of Re(z)
    numEntries = indicesB.size();
    indicesBReal = new int[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        indicesBReal[k] = 2 * indicesB[k];
    valuesBImag = new double[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        valuesBImag[k] = std::imag ( valuesB[k] );

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( imagRow, numEntries, valuesBImag, indicesBReal, firstTime ) );

    delete[] valuesBImag;
    valuesBImag = NULL;
    delete[] indicesBReal;
    indicesBReal = NULL;
    // -------------------------------------------------------------------
    // insert the coefficients Re(A) of Im(z)
    numEntries = indicesA.size();
    indicesAImag = new int[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        indicesAImag[k] = 2 * indicesA[k] + 1;
    valuesAReal = new double[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        valuesAReal[k] = std::real ( valuesA[k] );

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( imagRow, numEntries, valuesAReal, indicesAImag, firstTime ) );

    delete[] valuesAReal;
    valuesAReal = NULL;
    delete[] indicesAImag;
    indicesAImag = NULL;
    // -------------------------------------------------------------------
    // insert the coefficients -Re(B) of Im(z)
    numEntries = indicesB.size();
    indicesBImag = new int[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        indicesBImag[k] = 2 * indicesB[k] + 1;
    valuesBReal = new double[numEntries];
    for ( int k = 0; k < numEntries; k++ )
        valuesBReal[k] = -std::real ( valuesB[k] );

    TEUCHOS_ASSERT_INEQUALITY ( 0, <=, PutRow ( imagRow, numEntries, valuesBReal, indicesBImag, firstTime ) );

    delete[] valuesBReal;
    valuesBReal = NULL;
    delete[] indicesBImag;
    indicesBImag = NULL;
    // -------------------------------------------------------------------
}
// =============================================================================
int
Ginla::Komplex::
PutRow ( int Row,
         int & numIndices,
         double * values,
         int * indices,
         bool firstTime )
{
    if ( firstTime )
        return realMatrix_->InsertGlobalValues ( Row, numIndices, values, indices );
    else
        return realMatrix_->SumIntoGlobalValues ( Row, numIndices, values, indices );
}
// =============================================================================
Teuchos::RCP<Epetra_CrsMatrix>
Ginla::Komplex::getMatrix() const
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
Teuchos::RCP<Epetra_Vector>
Ginla::Komplex::
imagScalarProductCoeff ( const Teuchos::RCP<const ComplexVector> a ) const
{
    TEUCHOS_ASSERT ( a.is_valid_ptr() && !a.is_null() );
    TEUCHOS_ASSERT ( a->getMap()->isSameAs ( *ComplexMap_ ) );

    Teuchos::RCP<Epetra_Vector> compound = Teuchos::rcp ( new Epetra_Vector ( *RealMap_ ) );

    Teuchos::ArrayRCP<const double_complex> aView = a->get1dView();
    for ( unsigned int k=0; k<ComplexMap_->getNodeNumElements(); k++ )
    {
        ( *compound ) [2*k]   = -std::imag ( aView[k] );
        ( *compound ) [2*k+1] =  std::real ( aView[k] );
    }

    return compound;
}
// =============================================================================
