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

#ifndef VIO_IMAGE_WRITER_ABSTRACT_H
#define VIO_IMAGE_WRITER_ABSTRACT_H

#include "VIO_Typedefs.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Epetra_Vector.h>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>

namespace VIO {
namespace Image {
namespace Writer {

class Abstract
{
public:
    //! Default constructor.
    Abstract ( const std::string & filePath );

    //! Destructor
    virtual ~Abstract();

    virtual void
    write () const = 0; // pure virtual

    void
    setImageData ( const Epetra_MultiVector             & x,
                   const Teuchos::Tuple<unsigned int,2> & Nx,
                   const Point                          & h,
                   const Teuchos::Array<int>            & p = Teuchos::Array<int>(),
                   const Teuchos::Array<std::string>    & scalarsNames  = Teuchos::Array<std::string>()
                 );

    void
    setImageData ( const DoubleMultiVector              & x,
                   const Teuchos::Tuple<unsigned int,2> & Nx,
                   const Point                          & h,
                   const Teuchos::Array<int>            & p = Teuchos::Array<int>(),
                   const Teuchos::Array<std::string>    & scalarsNames = Teuchos::Array<std::string>()
                 );

    void
    setImageData ( const ComplexMultiVector             & z,
                   const Teuchos::Tuple<unsigned int,2> & Nx,
                   const Point                          & h,
                   const Teuchos::Array<int>            & p = Teuchos::Array<int>(),
                   const Teuchos::Array<std::string>    & scalarsNames = Teuchos::Array<std::string>()
                 );

    //! Add a parameter list to be stored in the field data section of the file.
    void
    addParameterList ( const Teuchos::ParameterList & problemParams );

    //! Add extra field data to be stored in the file.
    void
    addFieldData ( const Teuchos::Array<int> & array,
                   const std::string         & name );

protected:
    const std::string filePath_;
    const vtkSmartPointer<vtkImageData>   imageData_;

private:
};

  } // namespace Writer
 }
} // namespace VIO

#endif // VIO_IMAGE_WRITER_ABSTRACT_H
