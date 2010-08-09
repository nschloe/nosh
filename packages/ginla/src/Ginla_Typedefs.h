#ifndef GINLA_TYPEDEFS_H
#define GINLA_TYPEDEFS_H

#include <complex>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal


typedef Tpetra::Vector<double,Thyra::Ordinal> RealVector;
typedef Teuchos::Tuple<double,3> Point;

typedef std::complex<double> double_complex;
typedef Tpetra::Map<Thyra::Ordinal> ComplexMap;
typedef Tpetra::Vector<double_complex,Thyra::Ordinal> ComplexVector;
typedef Tpetra::MultiVector<double_complex,Thyra::Ordinal> ComplexMultiVector;
typedef Tpetra::CrsMatrix<double_complex,Thyra::Ordinal> ComplexMatrix;

const double_complex IM = double_complex( 0.0, 1.0 );

#endif // GINLA_TYPEDEFS_H