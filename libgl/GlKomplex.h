/*
 * GlKomplex.h
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl�mer
 */

#ifndef GLKOMPLEX_H_
#define GLKOMPLEX_H_

#include <Teuchos_RCP.hpp>
#include <Tpetra_Vector.hpp>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

typedef std::complex<double> double_complex;
typedef Tpetra::Vector<double_complex,Thyra::Ordinal> ComplexVector;


class GlKomplex
{

public:

  // Default constructor
  GlKomplex( const Teuchos::RCP<const Epetra_Comm>                  eComm,
	         const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > ComplexMap );

  // Destructor
  virtual
  ~GlKomplex();

  Teuchos::RCP<Epetra_Map>
  getRealMap() const;

  Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> >
  getComplexMap() const;

  //! Converts a real-valued vector to a complex-valued vector.
  Teuchos::RCP<ComplexVector>
  real2complex(const Epetra_Vector & x ) const;

  //! Converts a complex-valued vector to a real-valued vector.
  Teuchos::RCP<Epetra_Vector>
  complex2real( const ComplexVector &complexVec ) const;

  void
  initializeMatrix();

  void
  finalizeMatrix();

  void
  zeroOutMatrix();

  Teuchos::RCP<const Epetra_CrsMatrix>
  getMatrix() const;

  void
  updateRow( const int                            row,
             const std::vector<int>             & indicesA,
             const std::vector<double_complex>  & valuesA,
             const std::vector<int>             & indicesB,
             const std::vector<double_complex>  & valuesB,
             const bool                           firstTime
           );

private:

  Teuchos::RCP<Epetra_Map>
  createRealMap( const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > & ComplexMap ) const;

  Teuchos::RCP<const Teuchos::Comm<int> >
  create_CommInt( const Teuchos::RCP<const Epetra_Comm> &epetraComm );

  int
  PutRow( int Row, int & numIndices, double * values, int * indices, bool firstTime );

private:

  //! Communicators.
  const Teuchos::RCP<const Epetra_Comm>         EComm_;
  const Teuchos::RCP<const Teuchos::Comm<int> > TComm_;

  //! Maps for real- and complex-valued vectors.
  Teuchos::RCP<Epetra_Map>                         RealMap_;
  Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > ComplexMap_;

  Teuchos::RCP<Epetra_CrsMatrix> realMatrix_;

};

#endif /* GLKOMPLEX_H_ */
