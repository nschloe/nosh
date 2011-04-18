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
read( const Epetra_Comm                       & comm,
      Teuchos::RCP<Epetra_Vector>             & psi,
      Teuchos::RCP<Epetra_MultiVector>        & mvp,
      Teuchos::RCP<Epetra_Vector>             & thickness,
      Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
      Teuchos::ParameterList                  & parameterList
    )
{
  // Take two different fields with one component
  // instead of one field with two components. This works around
  // Ioss's inability to properly read psi_R, psi_Z as a complex variable.
  // (It can handle data_X, data_Y, data_Z though.)
  const unsigned int neq = 1;

  Teuchos::RCP<stk::mesh::fem::FEMMetaData> metaData =
      Teuchos::rcp( new stk::mesh::fem::FEMMetaData() );

  int numDim = 3;
  if (! metaData->is_FEM_initialized())
      metaData->FEM_initialize(numDim);

//   // attach fem data
//   size_t spatial_dimension = 3;
//   stk::mesh::DefaultFEM fem( *metaData, spatial_dimension );

  unsigned int field_data_chunk_size = 1001;
  Teuchos::RCP<stk::mesh::BulkData> bulkData =
      Teuchos::rcp( new stk::mesh::BulkData( stk::mesh::fem::FEMMetaData::get_meta_data(*metaData),
                                             MPI_COMM_WORLD,
                                             field_data_chunk_size
                                           )
                  );

  Teuchos::RCP<VectorFieldType> coordinatesField =
       Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "coordinates" ) );
  stk::io::set_field_role(*coordinatesField, Ioss::Field::ATTRIBUTE);

  // real part
  Teuchos::RCP<VectorFieldType> psir_field =
      Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "psi_R" ) );
  stk::mesh::put_field( *psir_field , metaData->node_rank() , metaData->universal_part(), neq );
  stk::io::set_field_role(*psir_field, Ioss::Field::TRANSIENT);

  // imaginary part
  Teuchos::RCP<VectorFieldType> psii_field =
      Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "psi_Z" ) );
  stk::mesh::put_field( *psii_field , metaData->node_rank() , metaData->universal_part(), neq );
  stk::io::set_field_role(*psii_field, Ioss::Field::TRANSIENT);

  // magnetic vector potential
  Teuchos::RCP<VectorFieldType> mvpXField =
       Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "AX" ) );
  stk::mesh::put_field( *mvpXField , metaData->node_rank() , metaData->universal_part(), neq );
  stk::io::set_field_role(*mvpXField, Ioss::Field::ATTRIBUTE);

  Teuchos::RCP<VectorFieldType> mvpYField =
       Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "AY" ) );
  stk::mesh::put_field( *mvpYField , metaData->node_rank() , metaData->universal_part(), neq );
  stk::io::set_field_role(*mvpYField, Ioss::Field::ATTRIBUTE);

  Teuchos::RCP<VectorFieldType> mvpZField =
       Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "AZ" ) );
  stk::mesh::put_field( *mvpZField , metaData->node_rank() , metaData->universal_part(), neq );
  stk::io::set_field_role(*mvpZField, Ioss::Field::ATTRIBUTE);

  Teuchos::RCP<VectorFieldType> thicknessField =
       Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "thickness" ) );
  stk::mesh::put_field( *thicknessField , metaData->node_rank() , metaData->universal_part(), neq );
  stk::io::set_field_role(*thicknessField, Ioss::Field::ATTRIBUTE);

  Teuchos::RCP<stk::io::util::MeshData> meshData =
      Teuchos::rcp( new stk::io::util::MeshData() );

  Ioss::Init::Initializer io;

  stk::io::util::create_input_mesh( "exodusii",
                                    fileName_,
                                    "",
                                    MPI_COMM_WORLD,
                                    *metaData,
                                    *meshData,
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
//   int index = -1; // Default to no restart
  int index = 1; // restart from the first step
  if ( comm.MyPID() == 0 )
      if ( index<1 )
        std::cout << "Restart Index not set. Not reading solution from exodus (" << index << ")"<< endl;
      else
        std::cout << "Restart Index set, reading solution time step: " << index << endl;

  stk::io::util::populate_bulk_data( *bulkData,
                                     *meshData,
                                     "exodusii",
                                     index
                                   );

  // add parameter
  meshData->m_region->field_add( Ioss::Field( "mu",
                                              Ioss::Field::REAL,
                                              "scalar",
                                              Ioss::Field::REDUCTION,
                                              1
                                            )
                               );

  bulkData->modification_end();

  //   coordinatesField = Teuchos::rcpFromRef( metaData->get_field<VectorFieldType>( std::string("coordinates") ) );

    // create the mesh with these specifications
    mesh = Teuchos::rcp( new Ginla::EpetraFVM::StkMesh( comm, metaData, bulkData, coordinatesField ) );

    // create the state
    psi       = this->createPsi_( mesh, psir_field, psii_field );
    mvp       = this->createMvp_( mesh, mvpXField, mvpYField, mvpZField );
    thickness = this->createThickness_( mesh, thicknessField );

    // These are vain attempts to find out whether thicknessField is actually empty.
//     const stk::mesh::FieldBase::RestrictionVector & restrictions = thicknessField->restrictions();
//     TEUCHOS_ASSERT( !restrictions.empty() );
//     std::cout << "max_size " << thicknessField->max_size(metaData->node_rank()) << std::endl;

    // Check of the thickness data is of any value. If not: ditch it.
    double norminf;
    thickness->NormInf( &norminf );
    if ( norminf < 1.0e-15 )
        thickness->PutScalar( 1.0 );

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
            const Teuchos::RCP<VectorFieldType>                 & psir_field,
            const Teuchos::RCP<VectorFieldType>                 & psii_field
          ) const
{
    // Get owned nodes.
    const std::vector<stk::mesh::Entity*> & ownedNodes = mesh->getOwnedNodes();

    // Create vector with this respective map.
    Teuchos::RCP<Epetra_Vector> psi = Teuchos::rcp( new Epetra_Vector( *mesh->getComplexMap() ) );

    // Fill the vector with data from the file
    int ind;
    for ( int k=0; k<ownedNodes.size(); k++ )
    {
        // real part
        double* psirVal = stk::mesh::field_data( *psir_field, *ownedNodes[k] );
        ind = 2*k;
        psi->ReplaceMyValues( 1, psirVal, &ind );

        // imaginary part
        double* psiiVal = stk::mesh::field_data( *psii_field, *ownedNodes[k] );
        ind = 2*k+1;
        psi->ReplaceMyValues( 1, psiiVal, &ind );
    }

    return psi;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::EpetraFVM::StkMeshReader::
createThickness_( const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
                  const Teuchos::RCP<VectorFieldType>                 & thickness_field
                ) const
{
    // Get owned nodes.
    const std::vector<stk::mesh::Entity*> & ownedNodes = mesh->getOwnedNodes();

    // Create vector with this respective map.
    Teuchos::RCP<Epetra_Vector> thickness = Teuchos::rcp( new Epetra_Vector( *mesh->getNodesMap() ) );

    TEUCHOS_ASSERT( !thickness_field.is_null() );

    // Fill the vector with data from the file
    for ( int k=0; k<ownedNodes.size(); k++ )
    {
        double* thicknessVal = stk::mesh::field_data( *thickness_field, *ownedNodes[k] );
        // Check if the field is actually there.
        if (thicknessVal == NULL)
        {
            std::cerr << "WARNING: Thickness value for node " << k << " not found.\n"
                      << "Probably there is no thickness field given with the state. Using default."
                      << std::endl;
            return Teuchos::null;
        }
        thickness->ReplaceMyValues( 1, thicknessVal, &k );
    }

    return thickness;
}
// =============================================================================
Teuchos::RCP<Epetra_MultiVector>
Ginla::EpetraFVM::StkMeshReader::
createMvp_( const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
            const Teuchos::RCP<const VectorFieldType>           & mvpXField,
            const Teuchos::RCP<const VectorFieldType>           & mvpYField,
            const Teuchos::RCP<const VectorFieldType>           & mvpZField
          ) const
{
    // Get owned nodes.
    const std::vector<stk::mesh::Entity*> & ownedNodes = mesh->getOwnedNodes();

    // Create vector with this respective map.
    int numComponents = 3;
    Teuchos::RCP<Epetra_MultiVector> mvp = Teuchos::rcp( new Epetra_MultiVector( *mesh->getNodesMap(), numComponents ) );

    TEUCHOS_ASSERT( !mvpXField.is_null() );
    TEUCHOS_ASSERT( !mvpYField.is_null() );
    TEUCHOS_ASSERT( !mvpZField.is_null() );

    // Fill the vector with data from the file
    for ( int k=0; k<ownedNodes.size(); k++ )
    {
        double* mvpXVal = stk::mesh::field_data( *mvpXField, *ownedNodes[k] );
        double* mvpYVal = stk::mesh::field_data( *mvpYField, *ownedNodes[k] );
        double* mvpZVal = stk::mesh::field_data( *mvpZField, *ownedNodes[k] );
        mvp->ReplaceMyValue( k, 0, *mvpXVal );
        mvp->ReplaceMyValue( k, 1, *mvpYVal );
        mvp->ReplaceMyValue( k, 2, *mvpZVal );
    }

    return mvp;
}
// =============================================================================
//
// Helper functions
//
void
Ginla::EpetraFVM::
StkMeshRead ( const Epetra_Comm                       & comm,
              const std::string                       & fileName,
              Teuchos::RCP<Epetra_Vector>             & psi,
              Teuchos::RCP<Epetra_MultiVector>        & mvp,
              Teuchos::RCP<Epetra_Vector>             & thickness,
              Teuchos::RCP<Ginla::EpetraFVM::StkMesh> & mesh,
              Teuchos::ParameterList                  & parameterList
            )
{
    StkMeshReader reader( fileName );

    reader.read( comm, psi, mvp, thickness, mesh, parameterList );

    return;
}
// =============================================================================
