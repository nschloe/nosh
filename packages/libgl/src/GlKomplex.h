/*
 * GlKomplex.h
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl\"omer
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
    GlKomplex ( const Teuchos::RCP<const Epetra_Comm>                  eComm,
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
    real2complex ( const Epetra_Vector & x ) const;

    //! Converts a complex-valued vector to a real-valued vector.
    Teuchos::RCP<Epetra_Vector>
    complex2real ( const ComplexVector &complexVec ) const;
    Teuchos::RCP<Epetra_Vector>
    complex2real ( const Teuchos::RCP<const ComplexVector> & complexVecPtr ) const;

    void
    finalizeMatrix();

    void
    zeroOutMatrix();

    //! Get read/write access to the matrix from outside.
    Teuchos::RCP<Epetra_CrsMatrix>
    getMatrix() const;



    /** \param row       row of the matrix that gets updated
      * \param indicesA  column indices of matrix \f$A\f$
      * \param valuesA   values in \f$A\f$ at \c indicesA
      * \param indicesB  column indices of matrix \f$B\f$
      * \param valuesA   values in \f$B\f$ at \c indicesB
      * \param firstTime whether or not the method gets called for the first time
      *
      * Of an equation system
      * \f[
      * A\psi + B \psi^* = b
      * \f]
      * where \f$A,B\in\mathbb{C}^{n\times n}\f$, \f$\psi, b\in\mathbb{C}^{n}\f$,
      * this routine constructs the corresponding real-valued equation system
      * \f[
      * \begin{pmatrix}
      * \Re{A}+\Re{B} & -\Im{A}+\Im{B}\\
      * \Im{A}+\Im{B} &  \Re{A}-\Re{B}
      * \end{pmatrix}
      * \begin{pmatrix}
      * \Re{\psi}\\
      * \Im{\psi}
      * \end{pmatrix}
      * =
      * \begin{pmatrix}
      * \Re{b}\\
      * \Im{b}
      * \end{pmatrix}
      * \f].
      */
    void
    updateRow ( const unsigned int                      row,
                const Teuchos::Array<int>             & indicesA,
                const Teuchos::Array<double_complex>  & valuesA,
                const Teuchos::Array<int>             & indicesB,
                const Teuchos::Array<double_complex>  & valuesB,
                const bool                              firstTime
              );

    /** Considering the term \f$a \psi + b \psi^*\f$ with
      * \f$a,b,\psi\in\mathbb{C}\f$, one has
      * \f[\begin{split}
      * \Re( a \psi + b \psi^* )
      * &= (\Re(a)+\Re(b))\Re(\psi) + (-\Im(a)+ \Im(b))\Im(\psi)\\
      * &= (\Re(a)+\Re(b),-\Im(a)+ \Im(b)) (\Re(\psi),\Im\psi)^{\mathrm{T}}
      * \end{split}
      * \f]
      * For given \f$a,b\f$, this function returns the vector \f$v\f$.
    */
    Teuchos::RCP<Epetra_Vector>
    getRealCoefficients ( const Teuchos::RCP<const ComplexVector> a,
                          const Teuchos::RCP<const ComplexVector> b
                        ) const;

    //! The term \f$\Im(\psi^*z)\f$ with \f$\psi,z\in\mathbb{C}^n\f$
    //! can be written as a scalar product \f$v^{\mathrm{T}} w\f$,
    //! where \f$w\in\mathbb{R}^{2n}\f$ is the real-valued representation
    //! of \f$z\f$.
    //! For given \f$\psi\f$, this function returns the corresponding
    //! \f$v\in\mathbb{R}^{2n}\f$.
    Teuchos::RCP<Epetra_Vector>
    imagScalarProductCoeff ( const Teuchos::RCP<const ComplexVector> a
                           ) const;

private:

    Teuchos::RCP<Epetra_Map>
    createRealMap ( const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > & ComplexMap ) const;

    Teuchos::RCP<const Teuchos::Comm<int> >
    create_CommInt ( const Teuchos::RCP<const Epetra_Comm> &epetraComm );

    int
    PutRow ( int Row, int & numIndices, double * values, int * indices, bool firstTime );

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
