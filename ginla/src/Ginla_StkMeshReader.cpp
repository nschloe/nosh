// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2010, 2011  Nico Schl\"omer
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
// =============================================================================
#include "Ginla_StkMeshReader.hpp"
#include "Ginla_StkMesh.hpp"

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

#include <Teuchos_VerboseObject.hpp>

#ifdef GINLA_TEUCHOS_TIME_MONITOR
  #include <Teuchos_TimeMonitor.hpp>
#endif

#include <stk_io/IossBridge.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <Ionit_Initializer.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Region.h>

namespace Ginla {
// =============================================================================
StkMeshReader::
StkMeshReader( const std::string & fileName ):
fileName_( fileName ),
#ifdef GINLA_TEUCHOS_TIME_MONITOR
readTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: StkMeshReader::read") ),
#endif
out_( Teuchos::VerboseObjectBase::getDefaultOStream() )
{
}
// =============================================================================
StkMeshReader::
~StkMeshReader()
{
}
// =============================================================================
void
StkMeshReader::
read( const Epetra_Comm      & comm,
      Teuchos::ParameterList & data
    )
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm(*readTime_);
#endif

  // Take two different fields with one component
  // instead of one field with two components. This works around
  // Ioss's inability to properly read psi_R, psi_Z as a complex variable.
  // (It can handle data_X, data_Y, data_Z though.)
  const unsigned int neq = 1;

  Teuchos::RCP<stk::mesh::fem::FEMMetaData> metaData =
      Teuchos::rcp( new stk::mesh::fem::FEMMetaData() );

  int numDim = 3;
  if ( !metaData->is_FEM_initialized() )
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
       Teuchos::rcpFromRef( metaData->declare_field<VectorFieldType>( "coordinates" ) );
  stk::io::set_field_role(*coordinatesField, Ioss::Field::ATTRIBUTE);

  // real part
  Teuchos::RCP<ScalarFieldType> psir_field =
      Teuchos::rcpFromRef( metaData->declare_field<ScalarFieldType>( "psi_R" ) );
  stk::mesh::put_field( *psir_field , metaData->node_rank(), metaData->universal_part() );
  stk::io::set_field_role( *psir_field, Ioss::Field::TRANSIENT );

  // imaginary part
  Teuchos::RCP<ScalarFieldType> psii_field =
      Teuchos::rcpFromRef( metaData->declare_field<ScalarFieldType>( "psi_Z" ) );
  stk::mesh::put_field( *psii_field , metaData->node_rank(), metaData->universal_part() );
  stk::io::set_field_role( *psii_field, Ioss::Field::TRANSIENT );

  // Magnetic vector potential.
  // Declare those fields as TRANSIENT to make sure they are written out to the
  // exodus file. Note, however, that this creates a large data overhead as the
  // same data is written out in each step although the data don't change.
  //
  // On 05/11/2011 02:44 PM, Gregory Sjaardema wrote:
  // For now, the B and C fields will have to be declared as TRANSIENT fields
  // since they are nodal fields on the universal set which currently doesn't
  // support attributes.  One of the stories I was supposed to work on for
  // this sprint was increasing the exodus capabilities supported by Ioss
  // (used underneath stk_io) and attributes on nodeblocks is one of the
  // things supported by exodus, but not by Ioss.  However, I got bogged down
  // by some big debugging and didn't finish that story.  Hopefully, it will
  // be done in June at which time you could use attribute fields on the
  // universal set...
  Teuchos::RCP<VectorFieldType> mvpField =
       Teuchos::rcpFromRef( metaData->declare_field< VectorFieldType >( "A" ) );
  stk::mesh::put_field( *mvpField, metaData->node_rank(), metaData->universal_part() );
  stk::io::set_field_role( *mvpField, Ioss::Field::ATTRIBUTE );

  // Thickness fields. Same as above.
  Teuchos::RCP<ScalarFieldType> thicknessField =
       Teuchos::rcpFromRef( metaData->declare_field< ScalarFieldType >( "thickness" ) );
  stk::mesh::put_field( *thicknessField, metaData->node_rank(), metaData->universal_part() );
  stk::io::set_field_role( *thicknessField, Ioss::Field::ATTRIBUTE );

  // initialize database communication
  Ioss::Init::Initializer io;

  // ---------------------------------------------------------------------------
  // The database is manually initialized to get explicit access to it
  // to be able to call set_field_separator(). This tells the database to
  // recognize 'AX', 'AY', 'AZ', as components of one vector field.
  // By default, the separator is '_' such that only fields 'A_*' would be
  // recognized as such.
  // 2011-09-26, Greg's mail:
  // "It is possible to manually open the database and create the Ioss::Region,
  // put that in the MeshData object and then call create_input_mesh(), but that
  // gets ugly.  I will try to come up with a better option... "
  std::string meshType = "exodusii";
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create( meshType,
                                                   fileName_,
                                                   Ioss::READ_MODEL,
                                                   MPI_COMM_WORLD
                                                 );
  TEST_FOR_EXCEPT_MSG( dbi == NULL || !dbi->ok(),
                               "ERROR: Could not open database '" << fileName_
                               << "' of type '" << meshType << "'."
                             );

  // set the vector field label separator
  dbi->set_field_separator(0);

  // create region to feed into meshData
  Ioss::Region *in_region = new Ioss::Region( dbi, "input_model" );
  // ---------------------------------------------------------------------------

  Teuchos::RCP<stk::io::MeshData> meshData =
      Teuchos::rcp( new stk::io::MeshData() );
  meshData->m_input_region = in_region;

  stk::io::create_input_mesh( meshType,
                              fileName_,
                              MPI_COMM_WORLD,
                              *metaData,
                              *meshData
                            );

  // define_input_fields() doesn't like the ATTRIBUTE fields; disable.
  // What was it good for anyways?
//  stk::io::define_input_fields( *meshData,
//                                *metaData
//                              );

//  stk::io::put_io_part_attribute( metaData->universal_part() );

  metaData->commit();

  stk::io::populate_bulk_data( *bulkData,
                               *meshData
                             );

//  // add parameter
//  meshData->m_region->field_add( Ioss::Field( "mu",
//                                              Ioss::Field::REAL,
//                                              "scalar",
//                                              Ioss::Field::REDUCTION,
//                                              1
//                                            )
//                               );
//
//  bulkData->modification_end();


  // Restart index to read solution from exodus file.
//   int index = -1; // Default to no restart
  int index = 1; // restart from the first step
//  if ( index<1 )
//      *out_ << "Restart Index not set. Not reading solution from exodus (" << index << ")"<< endl;
//  else
//      *out_ << "Restart Index set, reading solution time step: " << index << endl;

  stk::io::process_input_request( *meshData,
                                  *bulkData,
                                  index
                                );

  // create the mesh with these specifications
  Teuchos::RCP<Ginla::StkMesh> mesh = Teuchos::rcp( new Ginla::StkMesh( comm,
                                                      metaData,
                                                      bulkData,
                                                      coordinatesField
                                                    )
                     );

  data.setName( "data" );

  // add the mesh to the data list
  data.set( "mesh", mesh );

  // create the state
  data.set( "psi", this->createPsi_( mesh, psir_field, psii_field ) );
  data.set( "A", this->createMvp_( mesh, mvpField ) );

  // Check of the thickness data is of any value. If not: ditch it.
  Teuchos::RCP<Epetra_Vector> thickness = this->createThickness_( mesh, thicknessField );
  double norminf;
  thickness->NormInf( &norminf );
  if ( norminf < 1.0e-15 ) // assume that thickness wasn't present, fill with default value
      thickness->PutScalar( 1.0 );
  data.set( "thickness", thickness );

  // Add some dummy data.
  // TODO Replace by proper values.
  Teuchos::ParameterList problemParams;
  problemParams.setName( "Problem parameters" );
  problemParams.set( "mu", 0.0 );
  data.set( "Problem parameters", problemParams );

  // These are vain attempts to find out whether thicknessField is actually empty.
//     const stk::mesh::FieldBase::RestrictionVector & restrictions = thicknessField->restrictions();
//     TEUCHOS_ASSERT( !restrictions.empty() );
//     *out << "max_size " << thicknessField->max_size(metaData->node_rank()) << std::endl;

  return;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
StkMeshReader::
createPsi_( const Teuchos::RCP<const Ginla::StkMesh> & mesh,
            const Teuchos::RCP<ScalarFieldType>                 & psir_field,
            const Teuchos::RCP<ScalarFieldType>                 & psii_field
          ) const
{
    // Psi needs to have unique node IDs to be able to compute Norm2().
    // This is required in Belos.
    const std::vector<stk::mesh::Entity*> & ownedNodes = mesh->getOwnedNodes();

    // Create vector with this respective map.
    Teuchos::RCP<Epetra_Vector> psi = Teuchos::rcp( new Epetra_Vector( *mesh->getComplexNonOverlapMap() ) );

    // Fill the vector with data from the file.
    int ind;
    for ( unsigned int k=0; k<ownedNodes.size(); k++ )
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

    double r[1];
    psi->Norm1( r );
    TEST_FOR_EXCEPT_MSG( r[0]!=r[0] || r[0]>1.0e100,
                                 "The input data seems flawed. Abort." );

    return psi;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
StkMeshReader::
createThickness_( const Teuchos::RCP<const Ginla::StkMesh> & mesh,
                  const Teuchos::RCP<ScalarFieldType>                 & thickness_field
                ) const
{
    // Get overlap nodes.
    const std::vector<stk::mesh::Entity*> & overlapNodes = mesh->getOverlapNodes();

    // Create vector with this respective map.
    Teuchos::RCP<Epetra_Vector> thickness = Teuchos::rcp( new Epetra_Vector( *mesh->getNodesOverlapMap() ) );

    TEUCHOS_ASSERT( !thickness_field.is_null() );

    // Fill the vector with data from the file.
    for ( unsigned int k=0; k<overlapNodes.size(); k++ )
    {
        double* thicknessVal = stk::mesh::field_data( *thickness_field, *overlapNodes[k] );
        // Check if the field is actually there.
        if (thicknessVal == NULL)
        {
            *out_ << "WARNING: Thickness value for node " << k << " not found.\n"
                  << "Probably there is no thickness field given with the state. Using default."
                  << std::endl;
            return Teuchos::null;
        }
        int kk = int(k);
        thickness->ReplaceMyValues( 1, thicknessVal, &kk );
    }

    double r[1];
    thickness->Norm1( r );
    TEST_FOR_EXCEPT_MSG( r[0]!=r[0] || r[0]>1.0e100,
                                 "The input data seems flawed. Abort." );

    return thickness;
}
// =============================================================================
Teuchos::RCP<Epetra_MultiVector>
StkMeshReader::
createMvp_( const Teuchos::RCP<const Ginla::StkMesh> & mesh,
            const Teuchos::RCP<const VectorFieldType>           & mvpField
          ) const
{
    // Get overlap nodes.
    const std::vector<stk::mesh::Entity*> & overlapNodes = mesh->getOverlapNodes();

    // Create vector with this respective map.
    int numComponents = 3;
    Teuchos::RCP<Epetra_MultiVector> mvp = Teuchos::rcp( new Epetra_MultiVector( *mesh->getNodesOverlapMap(), numComponents ) );

    TEUCHOS_ASSERT( !mvpField.is_null() );

    // Fill the vector with data from the file.
    for ( unsigned int k=0; k<overlapNodes.size(); k++ )
    {
        double* mvpVal = stk::mesh::field_data( *mvpField, *overlapNodes[k] );
        // Check if the field is actually there.
        TEST_FOR_EXCEPT_MSG( mvpVal == NULL,
                                     "MVP value for node " << k << " not found.\n"
                                     << "Probably there is no MVP field given with the state."
                                   );

        mvp->ReplaceMyValue( k, 0, mvpVal[0] );
        mvp->ReplaceMyValue( k, 1, mvpVal[1] );
        mvp->ReplaceMyValue( k, 2, mvpVal[2] );
    }

    double r[3];
    mvp->Norm1( r );
    TEST_FOR_EXCEPT_MSG( r[0]!=r[0] || r[0]>1.0e100
                              || r[1]!=r[1] || r[1]>1.0e100
                              || r[2]!=r[2] || r[2]>1.0e100,
                                 "The input data seems flawed. Abort." );

    return mvp;
}
// =============================================================================
} // namespace Ginla
// =============================================================================
//
// Helper functions
//
void
Ginla::
StkMeshRead ( const Epetra_Comm      & comm,
              const std::string      & fileName,
              Teuchos::ParameterList & data
            )
{
    StkMeshReader reader( fileName );
    reader.read( comm, data );
    return;
}
// =============================================================================
