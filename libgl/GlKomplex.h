/*
 * GlKomplex.h
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schlšmer
 */

#ifndef GLKOMPLEX_H_
#define GLKOMPLEX_H_

#include <Teuchos_RCP.hpp>
#include <Tpetra_Vector.hpp>
#include <Epetra_Vector.h>
#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

typedef std::complex<double> double_complex;
typedef Tpetra::Vector<double_complex,Thyra::Ordinal> ComplexVector;


class GlKomplex
{

public:

  // Default constructor
  GlKomplex( const Teuchos::RCP<const Epetra_Comm> eComm );

  // Destructor
  virtual
  ~GlKomplex();

  //! Converts a real-valued vector to a complex-valued vector.
  Teuchos::RCP<ComplexVector>
  real2complex(const Epetra_Vector & x ) const;

  //! Converts a complex-valued vector to a real-valued vector.
  Teuchos::RCP<Epetra_Vector>
  complex2real( const ComplexVector &complexVec );

private:

  void
  createRealMap();

  Teuchos::RCP<const Teuchos::Comm<int> >
  create_CommInt( const Teuchos::RCP<const Epetra_Comm> &epetraComm );

private:

  //! Communicators.
  const Teuchos::RCP<const Epetra_Comm>         EComm_;
        Teuchos::RCP<const Teuchos::Comm<int> > TComm_;

  //! Maps for real- and complex-valued vectors.
  Teuchos::RCP<Epetra_Map>                         RealMap_;
  Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > ComplexMap_;

};

#endif /* GLKOMPLEX_H_ */
