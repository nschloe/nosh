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

#ifndef GINLA_EPETRAFVM_STKMESHREADER_H
#define GINLA_EPETRAFVM_STKMESHREADER_H
// =============================================================================
// includes
#include "Ginla_EpetraFVM_StkMesh.h"

#include <string>

#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
// =============================================================================
// typedefs
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
// =============================================================================
// forward declarations
class Epetra_Vector;
// =============================================================================
namespace Ginla {
namespace EpetraFVM {

class StkMeshReader
{
public:
    StkMeshReader( const std::string & fileName );

    virtual
    ~StkMeshReader();

    void
    read( const Epetra_Comm                       & comm,
          Teuchos::RCP<Epetra_Vector>             & psi,
          Teuchos::RCP<Epetra_MultiVector>        & mvp,
          Teuchos::RCP<Epetra_Vector>             & thickness,
          Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
          Teuchos::ParameterList                  & parameterList
        );

protected:
private:
    const std::string fileName_;

private:
    Teuchos::RCP<Epetra_Vector>
    createPsi_( const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
                const Teuchos::RCP<ScalarFieldType>                 & psir_field,
                const Teuchos::RCP<ScalarFieldType>                 & psii_field
              ) const;

    Teuchos::RCP<Epetra_Vector>
    createThickness_( const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
                      const Teuchos::RCP<ScalarFieldType>                 & thickness_field
                    ) const;

    Teuchos::RCP<Epetra_MultiVector>
    createMvp_( const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
                const Teuchos::RCP<const VectorFieldType>           & mvpField
              ) const;
};
// -----------------------------------------------------------------------------
// helper function
void
StkMeshRead ( const Epetra_Comm                       & comm,
              const std::string                       & fileName,
              Teuchos::RCP<Epetra_Vector>             & psi,
              Teuchos::RCP<Epetra_MultiVector>        & mvp,
              Teuchos::RCP<Epetra_Vector>             & thickness,
              Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
              Teuchos::ParameterList                  & parameterList
            );
// -----------------------------------------------------------------------------
} // namespace EpetraFVM
} // namespace Ginla
// =============================================================================
#endif // GINLA_EPETRAFVM_STKMESHREADER_H
