/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

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

#ifndef ABSTRACTIMAGEREADER_H
#define ABSTRACTIMAGEREADER_H

#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include <Teuchos_ParameterList.hpp>
#include <Tpetra_Vector.hpp>
#include <Thyra_OperatorVectorTypes.hpp> // For Thyra::Ordinal

typedef Tpetra::Vector<double,Thyra::Ordinal> DoubleVector;
typedef Tpetra::MultiVector<double,Thyra::Ordinal> DoubleMultiVector;
typedef Tpetra::Vector<std::complex<double>,Thyra::Ordinal> ComplexVector;
typedef Tpetra::MultiVector<std::complex<double>,Thyra::Ordinal> ComplexMultiVector;

typedef Teuchos::Tuple<int,2>      IntTuple;
typedef Teuchos::Tuple<unsigned,2> UIntTuple;
typedef Teuchos::Tuple<double,2>   DoubleTuple;

namespace VIO 
{
  namespace Reader 
  {

class Abstract
{
public:
    Abstract ( const std::string & filePath );

    virtual
    ~Abstract ();

    virtual void
    read ( Teuchos::RCP<ComplexMultiVector>              & z,
           Teuchos::Array<int>                           & p,
           UIntTuple                                     & dims,
           DoubleTuple                                   & origin,
           DoubleTuple                                   & spacing,
           Teuchos::ParameterList                        & fieldData,
           const Teuchos::RCP<const Teuchos::Comm<int> > & TComm
         ) const = 0; // purely virtual

protected:
    int
    extractIntValue ( const vtkSmartPointer<vtkDataArray> & array
                    ) const;

    Teuchos::Array<int>
    extractIntArray ( const vtkSmartPointer<vtkDataArray> & array
                    ) const;

    double
    extractDoubleValue ( const vtkSmartPointer<vtkDataArray> & array
                       ) const;

    Teuchos::Array<double>
    extractDoubleArray ( const vtkSmartPointer<vtkDataArray> & array
                       ) const;

protected:
    void
    processImageData ( const vtkSmartPointer<vtkImageData>           & imageData,
                       Teuchos::RCP<ComplexMultiVector>              & z,
                       Teuchos::Array<int>                           & p,
                       UIntTuple                                     & dims,
                       DoubleTuple                                   & origin,
                       DoubleTuple                                   & spacing,
                       Teuchos::ParameterList                        & fieldData,
                       const Teuchos::RCP<const Teuchos::Comm<int> > & TComm
                     ) const;

    Teuchos::ParameterList
    readFieldData ( const vtkSmartPointer<vtkImageData>  & imageData
                  ) const;

    std::string filePath_;

private:
};

  } // namespace Reader
} // namespace VIO

#endif // ABSTRACTIMAGEREADER_H
