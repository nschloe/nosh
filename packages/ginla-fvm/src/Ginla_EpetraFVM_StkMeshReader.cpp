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
    // Number of equations per node.
    unsigned int neq = 1;

    // Exodus files can contain more than one state.
    // 'index' tells which state is to be read.
    int index = 1;

    // initialize data sets
    Teuchos::RCP<stk::mesh::MetaData> metaData =
        Teuchos::rcp( new stk::mesh::MetaData(stk::mesh::fem_entity_rank_names() ) );
    unsigned int field_data_chunk_size = 1001;
    Teuchos::RCP<stk::mesh::BulkData> bulkData =
        Teuchos::rcp( new stk::mesh::BulkData(*metaData, MPI_COMM_WORLD , field_data_chunk_size ) );

    // read the coordinates of the nodes
    Teuchos::RCP<VectorFieldType> coordinates_field =
        Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "coordinates" ) );
    stk::io::set_field_role(*coordinates_field, Ioss::Field::ATTRIBUTE); // doesn't change for different time steps

    // read the real part
    Teuchos::RCP<VectorFieldType> psiR_field =
        Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "psi_R" ) );
    stk::mesh::put_field( *psiR_field , stk::mesh::Node , metaData->universal_part() , neq );
    stk::io::set_field_role(*psiR_field, Ioss::Field::TRANSIENT); // changes with steps

    // read the imaginary part
    Teuchos::RCP<VectorFieldType> psiI_field =
        Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "psi_Z" ) );
    stk::mesh::put_field( *psiI_field , stk::mesh::Node , metaData->universal_part() , neq );
    stk::io::set_field_role(*psiI_field, Ioss::Field::TRANSIENT); // changes with steps

    Teuchos::RCP<stk::io::util::MeshData> mesh_data = Teuchos::rcp( new stk::io::util::MeshData() );

#ifdef HAVE_MPI
    const Epetra_MpiComm& mpicomm = Teuchos::dyn_cast<const Epetra_MpiComm>( comm );
    MPI_Comm mcomm = mpicomm.Comm();
#else
    int mcomm = 1;
#endif

    // read the mesh and meta data from file
    Ioss::Init::Initializer io;
    std::string workingDirectory = "";
    stk::io::util::create_input_mesh( "exodusii",
                                      fileName_,
                                      workingDirectory,
                                      mcomm,
                                      *metaData,
                                      *mesh_data,
                                      false
                                    );

    std::map<int, stk::mesh::Part*> partVec;    // Element blocks
    std::map<std::string, stk::mesh::Part*> nsPartVec;  // Node Sets

    partVec[0] = & metaData->universal_part();
    stk::io::put_io_part_attribute(*partVec[0]);

    metaData->commit();

    // actual data is read here
    stk::io::util::populate_bulk_data(*bulkData, *mesh_data, "exodusii", index);
    bulkData->modification_end();

    // Create the mesh.
    mesh = Teuchos::rcp( new Ginla::EpetraFVM::StkMesh( comm, metaData, bulkData, coordinates_field ) );

    // create the state
    psi = this->createPsi_( comm, metaData, bulkData, psiR_field, psiI_field );

    // Add some dummy data.
    // TODO Replace by proper values.
    parameterList.set<double>( "mu", 0.0 );
    parameterList.set<double>( "scaling", 1.0 );

    return;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::EpetraFVM::StkMeshReader::
createPsi_( const Epetra_Comm & comm,
            const Teuchos::RCP<const stk::mesh::MetaData> & metaData,
            const Teuchos::RCP<const stk::mesh::BulkData> & bulkData,
            const Teuchos::RCP<VectorFieldType> & psiR_field,
            const Teuchos::RCP<VectorFieldType> & psiI_field
          ) const
{
    // Get owned nodes.
    std::vector<stk::mesh::Entity*> ownedNodes;
    stk::mesh::Selector select_owned_in_part = stk::mesh::Selector( metaData->universal_part() )
                                             & stk::mesh::Selector( metaData->locally_owned_part() );
    stk::mesh::get_selected_entities( select_owned_in_part,
                                      bulkData->buckets( stk::mesh::Node ),
                                      ownedNodes
                                    );

    // Create a map for real/imaginary out of this.
    int numDof = 2 * ownedNodes.size();
    std::vector<int> indices(numDof);
    for (int k=0; k < ownedNodes.size(); k++)
    {
        int globalNodeId = ownedNodes[k]->identifier() - 1;
        indices[2*k]   = 2*globalNodeId;
        indices[2*k+1] = 2*globalNodeId + 1;
    }
    Teuchos::RCP<Epetra_Map> complexMap =
        Teuchos::rcp(new Epetra_Map( -1, numDof, &(indices[0]), 0, comm) );

    // Create vector with this respective map.
    Teuchos::RCP<Epetra_Vector> psi = Teuchos::rcp( new Epetra_Vector( *complexMap ) );

    // Fill the vector with data from the file
    for ( int k=0; k<ownedNodes.size(); k++ )
    {
        double* sol = stk::mesh::field_data( *psiR_field, *ownedNodes[k] );
        int index = 2*k;
        psi->ReplaceMyValues( 1, sol, &index );

        sol = stk::mesh::field_data( *psiI_field, *ownedNodes[k] );
        index = 2*k+1;
        psi->ReplaceMyValues( 1, sol, &index );
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