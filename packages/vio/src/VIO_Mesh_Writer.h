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

#ifndef VIO_MESH_WRITER_H
#define VIO_MESH_WRITER_H
// =============================================================================
#include "VIO_Typedefs.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Epetra_Vector.h>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
// =============================================================================
namespace VIO {
  namespace Mesh {
    class Mesh;
  }
}
// =============================================================================
namespace VIO {
namespace Mesh {
class Writer
{
public:
    Writer ( const std::string & filePath );

    virtual
    ~Writer ();

    virtual void
    write () const; // pure virtual

    void
    setMesh( const VIO::Mesh::Mesh & mesh );

    void
    setValues( const Epetra_MultiVector          & x,
               const Teuchos::Array<std::string> & scalarsNames = Teuchos::Array<std::string>()
             );
               
    void
    setValues( const Tpetra::MultiVector<double> & x,
               const Teuchos::Array<std::string> & scalarsNames = Teuchos::Array<std::string>()
             );
             
    void
    setValues( const ComplexMultiVector          & z,
               const Teuchos::Array<std::string> & scalarsNames = Teuchos::Array<std::string>()
             );

//    //! Add a parameter list to be stored in the field data section of the file.
//     void
//     addParameterList ( const Teuchos::ParameterList & problemParams );
// 
//     //! Add extra field data to be stored in the file.
//     void
//     addFieldData ( const Teuchos::Array<int> & array,
//                    const std::string         & name );
  protected:
  private:
    const std::string filePath_;
    vtkSmartPointer<vtkUnstructuredGrid> vtkMesh_;
    Teuchos::RCP<VIO::Mesh::Mesh> mesh_;
};
}
}
// =============================================================================
#endif // VIO_MESH_WRITER_H