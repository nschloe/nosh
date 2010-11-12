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
#include "Ginla_EpetraFVM_StkMeshWriter.h"
#include "Ginla_EpetraFVM_StkMesh.h"

#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
// #include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/base/GetEntities.hpp>

// #include <stk_io/IossBridge.hpp>
#include <stk_io/util/UseCase_mesh.hpp>
// #include <Ionit_Initializer.h>
// =============================================================================
// typedefs
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
// =============================================================================
Ginla::EpetraFVM::StkMeshWriter::
StkMeshWriter( const std::string & fileName ):
fileName_( fileName ),
time_( 0 )
{
}
// =============================================================================
Ginla::EpetraFVM::StkMeshWriter::
~StkMeshWriter()
{
}
// =============================================================================
void
Ginla::EpetraFVM::StkMeshWriter::
write( const Epetra_Vector & psi,
       const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
       const Teuchos::ParameterList & parameterList
     )
{
    // create dummy mesh data
    Teuchos::RCP<stk::io::util::MeshData> meshData =
        Teuchos::rcp( new stk::io::util::MeshData() );

    // Get a native MPI comminicator object.
#ifdef HAVE_MPI
    const Epetra_MpiComm& mpicomm = Teuchos::dyn_cast<const Epetra_MpiComm&>( psi.Comm() );
    MPI_Comm mcomm = mpicomm.Comm();
#else
    int mcomm = 1;
#endif

    std::string meshExtension = "";
    std::string workingDirectory = "";
    // prepare the data for output
    stk::io::util::create_output_mesh( fileName_,
                                       meshExtension,
                                       workingDirectory,
                                       mcomm,
                                       *mesh->getBulkData(),
                                       *mesh->getMetaData(),
                                       *meshData
                                     );

    // Merge the state into the mesh.
    this->mergePsi_(  mesh->getMetaData(), mesh->getBulkData(), psi );

    // Write it.
    int out_step = stk::io::util::process_output_request( *meshData,
                                                          *mesh->getBulkData(),
                                                          time_++
                                                        );

    std::cout << "Ginla::EpetraFVM::StkMeshWriter::write:\n"
              << "\twriting time " << time << "\n"
              << "\tindex " << out_step << "\n"
              << "\tto file "<< fileName_
              << std::endl;

    return;
}
// =============================================================================
void
Ginla::EpetraFVM::StkMeshWriter::
mergePsi_( const Teuchos::RCP<stk::mesh::MetaData> & metaData,
           const Teuchos::RCP<stk::mesh::BulkData> & bulkData,
           const Epetra_Vector                     & psi
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

    VectorFieldType * psiR_field = metaData->get_field<VectorFieldType>( "psi_R" );
    VectorFieldType * psiI_field = metaData->get_field<VectorFieldType>( "psi_Z" );

    // Merge psi into the mesh.
    for (int k=0; k < ownedNodes.size(); k++)
    {
        // Extract real and imaginary part.
        double* psiR = stk::mesh::field_data( *psiR_field, *ownedNodes[k] );
        psiR[0] = psi[2*k];

        double* psiI = stk::mesh::field_data( *psiI_field, *ownedNodes[k] );
        psiI[0] = psi[2*k+1];
    }

    return;
}
// =============================================================================
//
// Helper functions
//
void
Ginla::EpetraFVM::
StkMeshWrite ( const std::string & fileName,
               const Epetra_Vector & psi,
               const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
               const Teuchos::ParameterList & parameterList
             )
{
    StkMeshWriter writer( fileName );
    writer.write( psi, mesh, parameterList );

    return;
}
// =============================================================================