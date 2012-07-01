// @HEADER
//
//    Reads stk_meshes.
//    Copyright (C) 2010--2012  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER
#ifndef CUANTICO_STKMESHREADER_H
#define CUANTICO_STKMESHREADER_H
// =============================================================================
// includes
#include "Cuantico_config.h"
#include "Cuantico_StkMesh.hpp"

#include <string>

#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>
#ifdef CUANTICO_TEUCHOS_TIME_MONITOR
  #include <Teuchos_Time.hpp>
#endif
#include <Teuchos_ParameterList.hpp>

// =============================================================================
// typedefs
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::mesh::Field<int> IntScalarFieldType ;
// =============================================================================
// forward declarations
class Epetra_Vector;
// =============================================================================
namespace Cuantico {

class StkMeshReader
{
public:
StkMeshReader( const std::string &fileName );

virtual
~StkMeshReader();

void
read( const Epetra_Comm &comm,
      Teuchos::ParameterList &data
      );

protected:
private:
const std::string fileName_;
const Teuchos::RCP<Teuchos::FancyOStream> out_;

private:

void
my_populate_bulk_data_(stk::mesh::BulkData &bulk_data,
                       Ioss::Region &region,
                       stk::mesh::fem::FEMMetaData &metaData);

Teuchos::RCP<Epetra_Vector>
complexfield2vector_( const Teuchos::RCP<const Cuantico::StkMesh> &mesh,
                      const Teuchos::RCP<ScalarFieldType> &realField,
                      const Teuchos::RCP<ScalarFieldType> &imagField
                      ) const;

Teuchos::RCP<Epetra_Vector>
scalarfield2vector_( const Teuchos::RCP<const Cuantico::StkMesh> &mesh,
                     const Teuchos::RCP<ScalarFieldType> &field
                     ) const;

Teuchos::RCP<Epetra_MultiVector>
createMvp_( const Teuchos::RCP<const Cuantico::StkMesh> &mesh,
            const Teuchos::RCP<const VectorFieldType> &mvpField
            ) const;

Teuchos::RCP<Epetra_MultiVector>
createMvpRZ_( const Teuchos::RCP<const Cuantico::StkMesh> &mesh,
              const Teuchos::RCP<const ScalarFieldType> &mvpFieldR,
              const Teuchos::RCP<const ScalarFieldType> &mvpFieldZ
              ) const;
};
// -----------------------------------------------------------------------------
// helper function
void
StkMeshRead( const Epetra_Comm &comm,
             const std::string &fileName,
             Teuchos::ParameterList &data
             );
// -----------------------------------------------------------------------------
} // namespace Cuantico
// =============================================================================
#endif // CUANTICO_STKMESHREADER_H
