/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl\"omer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef KOMPLEX2_DOUBLEMATRIX_H
#define KOMPLEX2_DOUBLEMATRIX_H

#include "Komplex2_Typedefs.h"

#include <Teuchos_RCP.hpp>

namespace Komplex2 {

class DoubleMatrix
{
  public:

  DoubleMatrix( const Teuchos::RCP<const ComplexMap> rangeMap,
                const Teuchos::RCP<const ComplexMap> domainMap
              );

  //! Getter for first matrix (corresponding to \f$\psi\f$).
  Teuchos::RCP<const ComplexMatrix>
  getMatrixA() const;

  //! Getter for second matrix (corresponding to \f$\\overline{psi}\f$).
  Teuchos::RCP<const ComplexMatrix>
  getMatrixB() const;

  void
  putALocalValues( unsigned int localRow,
                   const Teuchos::ArrayView<const Thyra::Ordinal> & columnIndices,
                   const Teuchos::ArrayView<const double_complex> & values
                 );

  void
  AsumIntoGlobalValues( unsigned int localRow,
                        const Teuchos::ArrayView<const Thyra::Ordinal> & columnIndices,
                        const Teuchos::ArrayView<const double_complex> & values
                      );

  void
  putBLocalValues( unsigned int localRow,
                   const Teuchos::ArrayView<const Thyra::Ordinal> & columnIndices,
                   const Teuchos::ArrayView<const double_complex> & values
                 );

  void
  BsumIntoGlobalValues( unsigned int localRow,
                        const Teuchos::ArrayView<const Thyra::Ordinal> & columnIndices,
                        const Teuchos::ArrayView<const double_complex> & values
                      );

  void
  finalize();

  protected:
  private:

  const Teuchos::RCP<ComplexMatrix> A_;
  const Teuchos::RCP<ComplexMatrix> B_;

};

}

#endif // KOMPLEX2_DOUBLEMATRIX_H
