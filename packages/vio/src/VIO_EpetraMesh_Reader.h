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

#ifndef VIO_EPETRAMESH_READER_H
#define VIO_EPETRAMESH_READER_H
// =============================================================================
#include "VIO_Reader_Abstract.h"
#include "VIO_Typedefs.h"

#include <string>
#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
// =============================================================================
// forward declarations
namespace VIO {
  namespace EpetraMesh {
      class Mesh;
  }
}
// =============================================================================
namespace VIO {
 namespace EpetraMesh {
class Reader:
   public VIO::Reader::Abstract
{
  public:
      Reader ( const std::string & filePath
             );

      ~Reader ();

      void
      read ( Teuchos::RCP<Epetra_Vector>           & z,
             Teuchos::RCP<EpetraMesh::Mesh>        & mesh,
             Teuchos::ParameterList                & fieldData,
             const Teuchos::RCP<const Epetra_Comm> & comm
           );

  protected:
  private:

    Teuchos::RCP<VIO::EpetraMesh::Mesh>
    extractMeshData_( const vtkSmartPointer<vtkUnstructuredGrid> & vtkMesh,
                      const Teuchos::RCP<const Epetra_Comm>      & comm
                    ) const;

    Teuchos::RCP<Epetra_Vector>
    extractStateData_ ( const vtkSmartPointer<vtkDataSet>     & vtkData,
                        const Teuchos::RCP<const Epetra_Comm> & comm
                      ) const;

    Teuchos::RCP<Epetra_Map>
    createComplexValuesMap_ ( const Epetra_Map  & nodesMap
                            ) const;
};

  // non-member function
  void
  read ( const Teuchos::RCP<const Epetra_Comm> & TComm,
         const std::string                     & filePath,
         Teuchos::RCP<Epetra_Vector>           & z,
         Teuchos::RCP<EpetraMesh::Mesh>        & mesh,
         Teuchos::ParameterList                & fieldData
        );
}
}
// =============================================================================
#endif // VIO_EPETRAMESH_READER_H
