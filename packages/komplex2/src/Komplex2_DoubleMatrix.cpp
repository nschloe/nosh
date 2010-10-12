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

#include "Komplex2_DoubleMatrix.h"

// ============================================================================
Komplex2::DoubleMatrix::
DoubleMatrix( const Teuchos::RCP<const ComplexMap> rangeMap,
              const Teuchos::RCP<const ComplexMap> domainMap
            ) :
            A_( Teuchos::rcp( new ComplexMatrix( rangeMap, domainMap, 0 ) ) ),
            B_( Teuchos::rcp( new ComplexMatrix( rangeMap, domainMap, 0 ) ) )
{
}
// ============================================================================            
Teuchos::RCP<const ComplexMatrix>
Komplex2::DoubleMatrix::
getMatrixA() const
{
  return A_;
}
// ============================================================================
Teuchos::RCP<const ComplexMatrix>
Komplex2::DoubleMatrix::
getMatrixB() const
{
  return B_;
}
// ============================================================================
void
Komplex2::DoubleMatrix::
putALocalValues( unsigned int localRow,
                 const Teuchos::ArrayView<const Thyra::Ordinal> & columnIndices,
                 const Teuchos::ArrayView<const double_complex> & values
               )
{
  if ( A_->isFillComplete() )  
     A_->replaceLocalValues( localRow, columnIndices, values );
  else
     A_->insertLocalValues( localRow, columnIndices, values );
}
// ============================================================================
void
Komplex2::DoubleMatrix::
AsumIntoGlobalValues( unsigned int localRow,
                      const Teuchos::ArrayView<const Thyra::Ordinal> & columnIndices,
                      const Teuchos::ArrayView<const double_complex> & values
                    )
{
  if ( A_->isFillComplete() )
     A_->sumIntoGlobalValues( localRow, columnIndices, values );
  else
     A_->insertLocalValues( localRow, columnIndices, values );
}
// ============================================================================
void
Komplex2::DoubleMatrix::
putBLocalValues( unsigned int localRow,
                 const Teuchos::ArrayView<const Thyra::Ordinal> & columnIndices,
                 const Teuchos::ArrayView<const double_complex> & values
               )
{
  if ( B_->isFillComplete() )  
     B_->replaceLocalValues( localRow, columnIndices, values );
  else
     B_->insertLocalValues( localRow, columnIndices, values );
}
// ============================================================================
void
Komplex2::DoubleMatrix::
BsumIntoGlobalValues( unsigned int localRow,
                     const Teuchos::ArrayView<const Thyra::Ordinal> & columnIndices,
                     const Teuchos::ArrayView<const double_complex> & values
                   )
{
  if ( B_->isFillComplete() )  
     B_->sumIntoGlobalValues( localRow, columnIndices, values );
  else
     B_->insertLocalValues( localRow, columnIndices, values );
}
// ============================================================================
void
Komplex2::DoubleMatrix::
finalize()
{ 
   // fill and optimize
   A_->fillComplete( Tpetra::DoOptimizeStorage );
   B_->fillComplete( Tpetra::DoOptimizeStorage );
}
// ============================================================================
