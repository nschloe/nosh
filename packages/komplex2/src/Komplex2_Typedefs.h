#ifndef KOMPLEX2_TYPEDEFS_H
#define KOMPLEX2_TYPEDEFS_H

#include <complex>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

typedef Thyra::Ordinal ORD;

typedef Tpetra::Vector<double,ORD> RealVector;
typedef Teuchos::Tuple<double,3> Point;

typedef std::complex<double> double_complex;
typedef Tpetra::Map<ORD> ComplexMap;
typedef Tpetra::Vector<double_complex,ORD> ComplexVector;
typedef Tpetra::MultiVector<double_complex,ORD> ComplexMultiVector;
typedef Tpetra::CrsMatrix<double_complex,ORD> ComplexMatrix;

const double_complex IM = double_complex( 0.0, 1.0 );

#endif // KOMPLEX2_TYPEDEFS_H
