/*
 * GlKomplex.h
 *
 *  Created on: Dec 16, 2009
 *      Author: Nico Schl\"omer
 */

#ifndef GINLA_KOMPLEX_LINEARPROBLEM_H_
#define GINLA_KOMPLEX_LINEARPROBLEM_H_

#include "Ginla_Typedefs.h"

// forward declarations
class Epetra_Comm;
class Epetra_Map;
class Epetra_CrsMatrix;
class Epetra_Vector;

namespace Ginla {
  namespace Komplex {
    class DoubleMatrix;
  }
}

namespace Ginla {

namespace Komplex {
  
class LinearProblem
{
public:

    // Default constructor
    LinearProblem ( const Teuchos::RCP<const Epetra_Comm>                  eComm,
                    const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > ComplexMap );

    // Destructor
    virtual
    ~LinearProblem();

    Teuchos::RCP<Epetra_Map>
    getRealMap() const;

    Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> >
    getComplexMap() const;

    //! Converts a real-valued vector to a complex-valued vector.
    Teuchos::RCP<ComplexVector>
    real2complex ( const Epetra_Vector & x ) const;

//     Teuchos::RCP<ComplexMultiVector>
//     real2complex ( const Epetra_MultiVector & x ) const;
    
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
    
    //! Essentially does the same thing as \c updateRow, just that the input argument
    //! is a whole \c DoubleMatrix.
    void
    update ( const Teuchos::RCP<const Ginla::Komplex::DoubleMatrix> AB,
             bool firstTime );

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
    updateRow ( const unsigned int                               row,
                const Teuchos::ArrayRCP<const Thyra::Ordinal>  & indicesA,
                const Teuchos::ArrayRCP<const double_complex>  & valuesA,
                const Teuchos::ArrayRCP<const Thyra::Ordinal>  & indicesB,
                const Teuchos::ArrayRCP<const double_complex>  & valuesB,
                const bool                                       firstTime
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
    //! For given \f$\psi\in\mathbb{C}^n\f$, this function returns the corresponding
    //! \f$v\in\mathbb{R}^{2n}\f$.
//     Teuchos::RCP<Epetra_Vector>
//     imagScalarProductCoeff ( const Teuchos::RCP<const ComplexVector> a
//                            ) const;

private:

    Teuchos::RCP<Epetra_Map>
    createRealMap_ ( const Teuchos::RCP<const Tpetra::Map<Thyra::Ordinal> > & ComplexMap ) const;

    Teuchos::RCP<const Teuchos::Comm<int> >
    create_CommInt_ ( const Teuchos::RCP<const Epetra_Comm> &epetraComm );

    int
    PutRow_ ( int      Row,
              int    & numIndices,
              double * values,
              int    * indices,
              bool     firstTime );

private:

    //! Communicators.
    const Teuchos::RCP<const Epetra_Comm>         EComm_;
    const Teuchos::RCP<const Teuchos::Comm<int> > TComm_;

    //! Maps for real- and complex-valued vectors.
    Teuchos::RCP<Epetra_Map>       RealMap_;
    Teuchos::RCP<const ComplexMap> ComplexMap_;

    Teuchos::RCP<Epetra_CrsMatrix> realMatrix_;

};

}

}
#endif /* GINLA_KOMPLEX_LINEARPROBLEM_H_ */
