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

#ifndef GINLA_EPETRAFVM_STKMESHWRITER_H
#define GINLA_EPETRAFVM_STKMESHWRITER_H
// =============================================================================
// includes
#include "Ginla_EpetraFVM_StkMesh.hpp"

#include <string>

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
// =============================================================================
namespace Ginla {
namespace EpetraFVM {

class StkMeshWriter
{
public:
    StkMeshWriter( const std::string & fileNameBase );

    virtual
    ~StkMeshWriter();

    void
    write( const Epetra_Vector & psi,
           const int             index,
           const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
           const Teuchos::ParameterList & parameterList
         );

protected:
private:
    const std::string fileNameBase_;
    int time_;
    const Teuchos::RCP<Teuchos::FancyOStream> out_;

private:
    void
    mergePsi_( const Teuchos::RCP<stk::mesh::fem::FEMMetaData> & metaData,
               const Teuchos::RCP<stk::mesh::BulkData>         & bulkData,
               const Epetra_Vector                             & psi
             ) const;

//    int
//    process_output_request_( stk::io::util::MeshData &mesh_data,
//                             stk::mesh::BulkData &bulk,
//                             double time,
//                             const Teuchos::ParameterList & parameterList,
//                             bool output_all_fields
//                           );
//
//    void
//    process_output_request_( Ioss::Region &region,
//                             stk::mesh::BulkData &bulk,
//                             int step,
//                             bool add_all_fields = false
//                           );

};
// -----------------------------------------------------------------------------
// helper function
void
StkMeshWrite ( const std::string & fileNameBase,
               const int           index,
               const Epetra_Vector & psi,
               const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
               const Teuchos::ParameterList & parameterList
            );
// -----------------------------------------------------------------------------
} // namespace EpetraFVM
} // namespace Ginla
// =============================================================================
#endif // GINLA_EPETRAFVM_STKMESHWRITER_H
