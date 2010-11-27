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
// =============================================================================
#include "Ginla_EpetraFVM_StkMeshReader.h"
#include "Ginla_EpetraFVM_StkMesh.h"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/util/UseCase_mesh.hpp>
#include <Ionit_Initializer.h>
// =============================================================================
// typedefs
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
// =============================================================================
Ginla::EpetraFVM::StkMeshReader::
StkMeshReader( const std::string & fileName ):
fileName_( fileName )
{
}
// =============================================================================
Ginla::EpetraFVM::StkMeshReader::
~StkMeshReader()
{
}
// =============================================================================
void
Ginla::EpetraFVM::StkMeshReader::
read( const Epetra_Comm & comm,
            Teuchos::RCP<Epetra_Vector> & psi,
            Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
            Teuchos::ParameterList & parameterList
    )
{
  const unsigned int neq = 2;

  Teuchos::RCP<stk::mesh::MetaData> metaData =
      Teuchos::rcp( new stk::mesh::MetaData( stk::mesh::fem_entity_rank_names() ) );
  unsigned int field_data_chunk_size = 1001;
  Teuchos::RCP<stk::mesh::BulkData> bulkData =
      Teuchos::rcp( new stk::mesh::BulkData( *metaData , MPI_COMM_WORLD , field_data_chunk_size ) );

  Teuchos::RCP<VectorFieldType> coordinatesField =
       Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "coordinates" ) );
  stk::io::set_field_role(*coordinatesField, Ioss::Field::ATTRIBUTE);

  Teuchos::RCP<VectorFieldType> solution_field =
      Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "psi" ) );
  stk::mesh::put_field( *solution_field , stk::mesh::Node , metaData->universal_part(), neq );
  stk::io::set_field_role(*solution_field, Ioss::Field::TRANSIENT);

  Teuchos::RCP<stk::io::util::MeshData> mesh_data =
      Teuchos::rcp( new stk::io::util::MeshData() );

  Ioss::Init::Initializer io;

  stk::io::util::create_input_mesh( "exodusii",
                                    fileName_,
                                    "",
                                    MPI_COMM_WORLD,
                                    *metaData,
                                    *mesh_data,
                                    false
                                  );

  stk::io::put_io_part_attribute( metaData->universal_part() );

//   // Set node sets
//   const stk::mesh::PartVector & all_parts = metaData->get_parts();
//   int eb = 0;
//   for (stk::mesh::PartVector::const_iterator i = all_parts.begin(); i != all_parts.end(); ++i)
//   {
//     stk::mesh::Part * const part = *i ;
// 
//     switch( part->primary_entity_rank() )
//     {
//       case stk::mesh::Element:
//       {
//           std::cout << "IOSS-STK: Element part found " << endl;
//           partVec[eb++] = part;
//           // Since Cubit likes to define numDim=3 always, use vertex
//           // count on top element block to figure out quad(tri) vs hex.
//           //   Needs to be fixed for Tets ro Shells
//           int numVerts = stk::mesh::get_cell_topology(*part)->vertex_count;
//           if (numVerts==4 || numVerts==3)
//               numDim=2;
//           else if (numVerts==8)
//               numDim=3;
//           else
//               TEST_FOR_EXCEPTION( true,
//                                   Teuchos::Exceptions::InvalidParameter,
//                                   std::endl << "Error!  IossSTKMeshStruct:  " <<
//                                   "Invalid vertex count from exodus mesh: " << numVerts << std::endl
//                                 );
//           std::cout << "IOSS-STK:  numDim =  " << numDim << endl;
//       }
//       break;
//       case stk::mesh::Node:
//           {
//             std:cout << "Mesh has Node Set ID: " << part->name() << endl;
//             nsPartVec[part->name()]=part;
//           }
//         break;
//       default:
//         break ;
//     }
//   }
// 
//   std::cout << "IOSS-STK: number of node sets = " << nsPartVec.size() << endl;

  metaData->commit();

  // Restart index to read solution from exodus file.
  int index = -1; // Default to no restart
  if ( index<1 )
    std::cout << "Restart Index not set. Not reading solution from exodus (" << index << ")"<< endl;
  else
    std::cout << "Restart Index set, reading solution time step: " << index << endl;

  stk::io::util::populate_bulk_data( *bulkData,
                                     *mesh_data,
                                     "exodusii",
                                     index
                                   );

  bulkData->modification_end();

  //   coordinatesField = Teuchos::rcpFromRef( metaData->get_field<VectorFieldType>( std::string("coordinates") ) );

    // create the mesh with these specifications
    mesh = Teuchos::rcp( new Ginla::EpetraFVM::StkMesh( comm, metaData, bulkData, coordinatesField ) );

    // create the state
    psi = this->createPsi_( mesh, solution_field );

    // Add some dummy data.
    // TODO Replace by proper values.
    parameterList.set<double>( "mu", 0.0 );
    parameterList.set<double>( "scaling", 1.0 );

    return;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::EpetraFVM::StkMeshReader::
createPsi_( const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
            const Teuchos::RCP<VectorFieldType>                 & psi_field
          ) const
{
    // Get owned nodes.
    const std::vector<stk::mesh::Entity*> & ownedNodes = mesh->getOwnedNodes();

    // Create vector with this respective map.
    Teuchos::RCP<Epetra_Vector> psi = Teuchos::rcp( new Epetra_Vector( *mesh->getComplexMap() ) );

    // Fill the vector with data from the file
    for ( int k=0; k<ownedNodes.size(); k++ )
    {
        double* psiVals = stk::mesh::field_data( *psi_field, *ownedNodes[k] );

        // TODO remove this hardcode
        psiVals[0] = 1.0;
        psiVals[1] = 0.0;

        int ind[2] = { 2*k, 2*k+1 };
        psi->ReplaceMyValues( 2, psiVals, ind );
    }

    return psi;
}
// =============================================================================
//
// Helper functions
//
void
Ginla::EpetraFVM::
StkMeshRead ( const Epetra_Comm & comm,
              const std::string & fileName,
              Teuchos::RCP<Epetra_Vector> & psi,
              Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
              Teuchos::ParameterList & parameterList
            )
{
    StkMeshReader reader( fileName );
    reader.read( comm, psi, mesh, parameterList );

    return;
}
// =============================================================================