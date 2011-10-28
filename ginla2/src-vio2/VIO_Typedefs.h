#ifndef VIO_TYPEDEFS_H
#define VIO_TYPEDEFS_H

#include <Tpetra_Vector.hpp>
#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

typedef Thyra::Ordinal ORD;

typedef std::complex<double> double_complex;
typedef Tpetra::Map<ORD> ComplexMap;

typedef Tpetra::Vector<double,ORD> DoubleVector;
typedef Tpetra::MultiVector<double,ORD> DoubleMultiVector;
typedef Tpetra::Vector<std::complex<double>,ORD> ComplexVector;
typedef Tpetra::MultiVector<std::complex<double>,ORD> ComplexMultiVector;

typedef Teuchos::Tuple<int,2>          IntTuple;
typedef Teuchos::Tuple<unsigned int,2> UIntTuple;
typedef Teuchos::Tuple<double,2>       DoubleTuple;

typedef Teuchos::Tuple<double,3> Point;

#endif // VIO_TYPEDEFS_H


