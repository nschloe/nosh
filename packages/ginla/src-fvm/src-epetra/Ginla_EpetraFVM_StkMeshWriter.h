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
#include "Ginla_EpetraFVM_StkMesh.h"

#include <string>

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
// =============================================================================
// forward declarations
class Epetra_Vector;
// =============================================================================
namespace Ginla {
namespace EpetraFVM {

class StkMeshWriter
{
public:
    StkMeshWriter( const std::string & fileName );

    virtual
    ~StkMeshWriter();

    void
    write( const Epetra_Vector & psi,
           const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
           const Teuchos::ParameterList & parameterList
         );

protected:
private:
    const std::string fileName_;
    int time_;

private:
    void
    mergePsi_( const Teuchos::RCP<stk::mesh::MetaData> & metaData,
               const Teuchos::RCP<stk::mesh::BulkData> & bulkData,
               const Epetra_Vector                     & psi
             ) const;
};
// -----------------------------------------------------------------------------
// helper function
void
StkMeshWrite ( const std::string & fileName,
               const Epetra_Vector & psi,
               const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
               const Teuchos::ParameterList & parameterList
            );
// -----------------------------------------------------------------------------
} // namespace EpetraFVM
} // namespace Ginla
// =============================================================================
#endif // GINLA_EPETRAFVM_STKMESHWRITER_H
