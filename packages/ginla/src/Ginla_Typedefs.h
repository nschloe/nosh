#ifndef GINLA_TYPEDEFS_H
#define GINLA_TYPEDEFS_H

#include <complex>
#include <Tpetra_Vector.hpp>
#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

typedef std::complex<double> double_complex;
typedef Tpetra::Vector<double_complex,Thyra::Ordinal> ComplexVector;
typedef Tpetra::MultiVector<double_complex,Thyra::Ordinal> ComplexMultiVector;

#endif // GINLA_TYPEDEFS_H