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
#include "Ginla_EpetraFVM_StkMeshWriter.hpp"
#include "Ginla_EpetraFVM_StkMesh.hpp"

#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_VerboseObject.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_io/IossBridge.hpp>

// =============================================================================
// typedefs
typedef stk::mesh::Field<double> ScalarFieldType;
// =============================================================================
namespace Ginla {
namespace EpetraFVM {
// =============================================================================
StkMeshWriter::
StkMeshWriter( const std::string & fileNameBase ):
fileNameBase_( fileNameBase ),
time_( 0 ),
out_( Teuchos::VerboseObjectBase::getDefaultOStream() )
{
}
// =============================================================================
StkMeshWriter::
~StkMeshWriter()
{
}
// =============================================================================
void
StkMeshWriter::
write( const Epetra_Vector                                 & psi,
       const int                                             step,
       const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
       const Teuchos::ParameterList                        & parameterList
     )
{

    // Get a native MPI comminicator object.
#ifdef HAVE_MPI
    const Epetra_MpiComm& mpicomm = Teuchos::dyn_cast<const Epetra_MpiComm>( psi.Map().Comm() );
    MPI_Comm mcomm = mpicomm.Comm();
#else
    int mcomm = 1;
#endif

    // Create file name.
    // For ParaView to read the files as a sequence, they must have this exact
    // format;
    // see <http://www.paraview.org/Wiki/Restarted_Simulation_Readers#Exodus>.

    std::stringstream meshExtension;
    if ( step==0 )
        meshExtension << "e";
    else
        meshExtension << "e-s." << std::setw(4) << std::setfill('0') << step;

    std::string workingDirectory = "";

    std::stringstream fileName;
    fileName << fileNameBase_ << "." << meshExtension.str();

    // prepare the data for output
    stk::io::create_output_mesh( fileName.str(),
                                 mcomm,
                                 *mesh->getBulkData(),
                                 *mesh->getMeshData()
                               );

    stk::io::define_output_fields( *mesh->getMeshData(),
                                   *mesh->getMetaData()
                                 );

    // Merge the state into the mesh.
//     mesh->getBulkData()->modification_begin();
    this->mergePsi_( mesh->getMetaData(), mesh->getBulkData(), psi );
//     mesh->getBulkData()->modification_end();

    // Write it.
    double time = step;
    int out_step = stk::io::process_output_request( *mesh->getMeshData(),
                                                    *mesh->getBulkData(),
                                                    time
                                                  );

    *out_ << "Ginla::EpetraFVM::StkMeshWriter::write:\n"
          << "\twriting time " << step << "\n"
          << "\tindex " << out_step << "\n"
          << "\tto file " //<< fileName.str()
          << std::endl;

    return;
}
// =============================================================================
void
StkMeshWriter::
mergePsi_( const Teuchos::RCP<stk::mesh::fem::FEMMetaData> & metaData,
           const Teuchos::RCP<stk::mesh::BulkData>         & bulkData,
           const Epetra_Vector                             & psi
         ) const
{
    // Get owned nodes.
    std::vector<stk::mesh::Entity*> ownedNodes;
    stk::mesh::Selector select_owned_in_part = stk::mesh::Selector( metaData->universal_part() )
                                             & stk::mesh::Selector( metaData->locally_owned_part() );
    stk::mesh::get_selected_entities( select_owned_in_part,
                                      bulkData->buckets( metaData->node_rank() ),
                                      ownedNodes
                                    );

    ScalarFieldType * psiR_field = metaData->get_field<ScalarFieldType>( "psi_R" );
    ScalarFieldType * psiI_field = metaData->get_field<ScalarFieldType>( "psi_Z" );

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
// This is a copy of $TRILINOS/packages/stk/stk_io/stk_io/util/UseCase_mesh.cpp,
// enhanced with parameter handling
//int
//Ginla::EpetraFVM::StkMeshWriter::
//process_output_request_( stk::io::util::MeshData &mesh_data,
//                         stk::mesh::BulkData &bulk,
//                         double time,
//                         const Teuchos::ParameterList & parameterList,
//                         bool output_all_fields
//                       )
//{
//   Ioss::Region &region = *(mesh_data.m_region);
//
//   region.begin_mode(Ioss::STATE_TRANSIENT);
//
//   int out_step = region.add_state(time);
//
////    // add the parameters
////    for ( Teuchos::ParameterList::ConstIterator item = parameterList.begin();
////          item != parameterList.end();
////          ++item )
////    {
////        std::string name =  parameterList.name(item);
////        std::vector<double>  val(1);
////        parameterList.entry(item).getValue( &val[0] );
////        *out << "Writing \"" << name << "\" with value \"" << val[0] << "\"." << std::endl;
////        region.put_field_data( name, val );
////    }
//
//  this->process_output_request_(region, bulk, out_step, output_all_fields);
//  region.end_mode(Ioss::STATE_TRANSIENT);
//
//  return out_step;
//}
//// =============================================================================
//// This is a copy of $TRILINOS/packages/stk/stk_io/stk_io/util/UseCase_mesh.cpp,
//// enhanced with parameter handling
//void
//Ginla::EpetraFVM::StkMeshWriter::
//process_output_request_( Ioss::Region &region,
//                         stk::mesh::BulkData &bulk,
//                         int step,
//                         bool output_all_fields
//                       )
//{
//
//  region.begin_state(step);
//  // Special processing for nodeblock (all nodes in model)...
//  const stk::mesh::MetaData & meta = bulk.mesh_meta_data();
//
//  stk::io::util::put_field_data( bulk,
//                                 meta.universal_part(),
//                                 stk::mesh::Node,
//                                 region.get_node_blocks()[0],
//                                 Ioss::Field::Field::TRANSIENT,
//                                 output_all_fields
//                               );
//
//  const stk::mesh::PartVector & all_parts = meta.get_parts();
//  for ( stk::mesh::PartVector::const_iterator
//          ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {
//
//    stk::mesh::Part * const part = *ip;
//
//    // Check whether this part should be output to results database.
//    if (stk::io::is_part_io_part(*part)) {
//      // Get Ioss::GroupingEntity corresponding to this part...
//      Ioss::GroupingEntity *entity = region.get_entity(part->name());
//      if (entity != NULL) {
//        if (entity->type() == Ioss::ELEMENTBLOCK) {
//          stk::io::util::put_field_data( bulk,
//                                         *part,
//                                         stk::mesh::fem_entity_rank( part->primary_entity_rank()),
//                                         entity,
//                                         Ioss::Field::Field::TRANSIENT,
//                                         output_all_fields
//                                       );
//        }
//      }
//    }
//  }
//  region.end_state(step);
//}
// =============================================================================
} // namespace EpetraFVM
} // namespace Ginla
//
// Helper functions
//
void
Ginla::EpetraFVM::
StkMeshWrite ( const std::string   & fileNameBase,
               const int             index,
               const Epetra_Vector & psi,
               const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
               const Teuchos::ParameterList & parameterList
             )
{
    StkMeshWriter writer( fileNameBase );
    writer.write( psi, index, mesh, parameterList );

    return;
}
// =============================================================================