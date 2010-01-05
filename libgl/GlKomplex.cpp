/*
 * GlKomplex.cpp
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl�mer
 */

#include "GlKomplex.h"

#include <Epetra_IntSerialDenseVector.h>
#include <Epetra_Map.h>

#include <Teuchos_DefaultComm.hpp> // for Teuchos::SerialComm

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif // HAVE_MPI

// abbreviate the complex type name
typedef std::complex<double> double_complex;

// =============================================================================
// Default constructor
GlKomplex::GlKomplex( const Teuchos::RCP<const Epetra_Comm> eComm ) :
  EComm_( eComm ),
  TComm_( create_CommInt(eComm) )
{
  // TODO There is (until now?) no way to convert a Teuchos::Comm (of psi)
  // to an Epetra_Comm (of the real valued representation of psi), so the
  // Epetra_Comm has to be generated explicitly, and two communicators are kept
  // side by side all the time. One must make sure that the two are actually
  // equivalent, which can be checked by Thyra's conversion method create_Comm.
  // TODO Is is actually necessary to have equivalent communicators on the
  // real-valued and the complex-valued side?
  // How to compare two communicators anyway?

  // create fitting Tpetra::Comm
  // TODO Move into initializer.
  // TODO Follow up on the development of Trilinos for templatized comms.
//  TComm_ = create_CommInt(EComm_);

//  // define maps
//  // psi->getMap() returns a CONST map
//  ComplexMap_ = Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> >(psi->getMap());
//
//  // get the map for the real values
//  makeRealMap(ComplexMap_);

//  // set the number of local elements
//  NumMyElements_ = RealMap_->NumMyElements();

}
// =============================================================================
GlKomplex::~GlKomplex()
{
}
// =============================================================================
// converts a real-valued vector to a complex-valued psi vector
Teuchos::RCP<ComplexVector>
GlKomplex::real2complex(const Epetra_Vector & x ) const
{
//  if( !x.Map().SameAs(RealMap_) ) {
//    // recreate maps
//    RealMap_    = x.Map();
//    ComplexMap_ = createComplexMap();
//  }

  TEST_FOR_EXCEPTION( !RealMap_.is_valid_ptr() || RealMap_.is_null(),
                      std::logic_error,
                      "RealMap_ not properly initialized." );

  TEST_FOR_EXCEPTION( !x.Map().SameAs(*RealMap_),
                      std::logic_error,
                      "Maps for real-valued vectors do not coincide. "
                      << "Check, for example, the number of elements "
                      << "(" << x.Map().NumGlobalElements() << " for x vs. "
                      << RealMap_->NumGlobalElements() << " for RealMap_).");

  TEST_FOR_EXCEPTION( !ComplexMap_.is_valid_ptr() || ComplexMap_.is_null(),
                      std::logic_error,
                      "ComplexMap_ has not been properly initialized." );

  Teuchos::RCP<ComplexVector> z = Teuchos::rcp( new ComplexVector(ComplexMap_) );

  // TODO: parallelize
  for (unsigned int k=0; k < z->getGlobalLength(); k++) {
      double_complex c = double_complex(x[2 * k], x[2 * k + 1]);
      z->replaceGlobalValue(k, c);
  }
  return z;
}
// =============================================================================
// converts a real-valued vector to a complex-valued psi vector
Teuchos::RCP<Epetra_Vector>
GlKomplex::complex2real( const ComplexVector &complexVec )
{
  if(    !ComplexMap_.is_valid_ptr()
      ||  ComplexMap_.is_null()
      || !complexVec.getMap()->isSameAs(*ComplexMap_) ) {
    // recreate maps
    ComplexMap_ = complexVec.getMap();
    createRealMap();
  }

  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp( new Epetra_Vector(*RealMap_) );
  // TODO: parallelize
  Teuchos::ArrayRCP<const double_complex> complexVecView = complexVec.get1dView();
  int numComplexUnknowns = ComplexMap_->getGlobalNumElements();

  for (int k = 0; k < numComplexUnknowns; k++) {
      x->ReplaceGlobalValue( 2*k  , 0, real(complexVecView[k]) );
      x->ReplaceGlobalValue( 2*k+1, 0, imag(complexVecView[k]) );
  }
  return x;
}
// =============================================================================
void
GlKomplex::createRealMap()
{
  TEST_FOR_EXCEPTION( !ComplexMap_.is_valid_ptr() || ComplexMap_.is_null(),
                      std::logic_error,
                      "ComplexMap_ has not been properly initialized." );

  // get view for the global indices of the global elements
  Teuchos::ArrayView<const Thyra::Ordinal> myComplexGIDs =
                                             ComplexMap_->getNodeElementList();

  // Construct the map in such a way that all complex entries on processor K
  // are split up into real and imaginary part, which will both reside on
  // processor K again.
  unsigned int numMyComplexElements = myComplexGIDs.size();
  unsigned int numMyRealElements    = 2*numMyComplexElements;
  Epetra_IntSerialDenseVector myRealGIDs( numMyRealElements );
  for (unsigned int i = 0; i < numMyComplexElements; i++) {
	  myRealGIDs[2*i  ] = 2 * myComplexGIDs[i];
	  myRealGIDs[2*i+1] = 2 * myComplexGIDs[i] + 1;
  }

  RealMap_ = Teuchos::rcp( new Epetra_Map(numMyRealElements,
		                                  myRealGIDs.Length(),
				                          myRealGIDs.Values(),
                                          ComplexMap_->getIndexBase(),
                                          *EComm_)
                         );
  return;
}
// =============================================================================
//void
//GlKomplex::createComplexMap()
//{
//  TEST_FOR_EXCEPTION( !RealMap_.is_valid_ptr() || RealMap_.is_null(),
//                      std::logic_error,
//                      "RealMap_ has not been properly initialized." );
//
//  n = RealMap_->NumGlobalElements();
//  TEST_FOR_EXCEPTION( !n%2,
//                      std::logical_error,
//                      "Number of elements in RealMap is odd (" << n
//                      << "). Can't create ComplexMap_.");
//
//  int numComplexGlobalElements =  n / 2;
//
//  Epetra_IntSerialDenseVector realMapGIDs(numRealGlobalElements);
//  Teuchos::ArrayView<const Thyra::Ordinal> myGlobalElements;
//  =
//                                             ComplexMap_->getNodeElementList();
//
//  // Construct the map in such a way that all complex entries on processor K
//  // are split up into real and imaginary part, which will both reside on
//  // processor K again.
//
//  for (unsigned int k=0; k<n; i++) {
//      myGlobalElements[k]
//      realMapGIDs[2*i  ] = 2 * myGlobalElements[i];
//      realMapGIDs[2*i+1] = 2 * myGlobalElements[i] + 1;
//  }
//
//  Complex_ = Teuchos::rcp( new Tpetra::Map(numRealGlobalElements,
//                                           elementList,
//                                           RealMap_->IndexBase(),
//                                           TComm_,
//
//                                           realMapGIDs.Length(),
//                                          realMapGIDs.Values(),
//                                          ComplexMap_->getIndexBase(),
//                                          *EComm_)
//                         );
//  return;
//}
// =============================================================================
Teuchos::RCP<const Teuchos::Comm<int> >
GlKomplex::create_CommInt( const Teuchos::RCP<const Epetra_Comm> &epetraComm )
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
  Teuchos::RCP<const Epetra_SerialComm>
    serialEpetraComm = rcp_dynamic_cast<const Epetra_SerialComm>(epetraComm);
  if( serialEpetraComm.get() ) {
    Teuchos::RCP<const Teuchos::SerialComm<int> >
      serialComm = rcp(new Teuchos::SerialComm<int>());
    set_extra_data( serialEpetraComm, "serialEpetraComm", Teuchos::inOutArg(serialComm) );
    return serialComm;
  }
#endif // HAVE_MPI

  // If you get here then the conversion failed!
  return Teuchos::null;
}
// =============================================================================
