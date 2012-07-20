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
// =============================================================================
#include "Nosh_StkMeshReader.hpp"
#include "Nosh_StkMesh.hpp"

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

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  #include <Teuchos_TimeMonitor.hpp>
#endif

// For the my_populate_bulk_data_.
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>

#include <stk_io/MeshReadWriteUtils.hpp>
#include <Ionit_Initializer.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Region.h>

#ifdef HAVE_MPI
// Rebalance
#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance_utils/RebalanceUtils.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>
#endif

namespace Nosh {
// =============================================================================
StkMeshReader::
StkMeshReader(const std::string &fileName,
              const int index) :
  fileName_( fileName ),
  index_( index ),
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
read( const Epetra_Comm &comm,
      Teuchos::ParameterList &data
      )
{
#ifdef HAVE_MPI
  const Epetra_MpiComm &mpicomm =
    Teuchos::dyn_cast<const Epetra_MpiComm>(comm);
  MPI_Comm mcomm = mpicomm.Comm();
#else
  int mcomm = 1;
#endif

  // Take two different fields with one component
  // instead of one field with two components. This works around
  // Ioss's inability to properly read psi_R, psi_Z as a complex variable.
  // (It can handle data_X, data_Y, data_Z though.)
  //const unsigned int neq = 1;

  Teuchos::RCP<stk::mesh::fem::FEMMetaData> metaData =
    Teuchos::rcp( new stk::mesh::fem::FEMMetaData() );

  int numDim = 3;
  if ( !metaData->is_FEM_initialized() )
    metaData->FEM_initialize( numDim );

//   // attach fem data
//   size_t spatial_dimension = 3;
//   stk::mesh::DefaultFEM fem( *metaData, spatial_dimension );

  // ---------------------------------------------------------------------------
  // initialize database communication
  Ioss::Init::Initializer io;

  // If the file is serial, read it with process 0 and embed it
  // in the multiproc context. Load balancing is done later anyways.
  MPI_Comm readerComm = mcomm;
#ifdef HAVE_MPI
  bool fileIsSerial = fileName_.substr(fileName_.find_last_of(".") + 1) == "e";
  if (fileIsSerial && comm.NumProc()>1)
  {
    // reader process
    int readerProc[1] = {0};

    // Get the group under mcomm.
    MPI_Group groupWorld;
    MPI_Comm_group(mcomm, &groupWorld);
    // Create the new group.
    MPI_Group peZero;
    MPI_Group_incl(groupWorld, 1, readerProc, &peZero);
    // Create the new communicator.
    MPI_Comm_create(mcomm, peZero, &readerComm);
  }
#endif

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
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(meshType,
                                                  fileName_,
                                                  Ioss::READ_MODEL,
                                                  readerComm);
  TEUCHOS_TEST_FOR_EXCEPT_MSG( dbi == NULL || !dbi->ok(),
                       "ERROR: Could not open database '" << fileName_
                       << "' of type '" << meshType << "'.");

  // set the vector field label separator
  dbi->set_field_separator( 0 );

  // create region to feed into meshData
  Ioss::Region *in_region = new Ioss::Region( dbi, "input_model" );
  // ---------------------------------------------------------------------------

  Teuchos::RCP<stk::io::MeshData> meshData =
    Teuchos::rcp( new stk::io::MeshData() );
  meshData->m_input_region = in_region;

  // This checks the existence of the file, checks to see if we can open it,
  // builds a handle to the region and puts it in mesh_data (in_region),
  // and reads the metaData into metaData.
  stk::io::create_input_mesh(meshType,
                             fileName_,
                             readerComm,
                             *metaData,
                             *meshData);
  // ---------------------------------------------------------------------------
  // Setup field data.
  // The work set size could probably be improved.
  // Check out what Albany does here.
  unsigned int field_data_chunk_size = 1001;
  Teuchos::RCP<stk::mesh::BulkData> bulkData =
    Teuchos::rcp(new stk::mesh::BulkData(stk::mesh::fem::FEMMetaData::
                                         get_meta_data( *metaData ),
                                         mcomm,
                                         field_data_chunk_size));

  // As of now (2012-03-21) there is no way to determine which fields are
  // actually present in the file. The only thing we can do is to declare them,
  // and check if they're full of zeros in the end.
  Teuchos::RCP<VectorFieldType> coordinatesField =
    Teuchos::rcpFromRef(metaData->declare_field<VectorFieldType>(
                           "coordinates"));
  stk::mesh::put_field(*coordinatesField,
                       metaData->node_rank(),
                       metaData->universal_part(),
                       numDim);
  stk::io::set_field_role(*coordinatesField,
                          Ioss::Field::MESH);

  //Teuchos::RCP<IntScalarFieldType> procRankField =
  //  Teuchos::rcpFromRef(metaData->declare_field<IntScalarFieldType>( "proc_rank" ));
  //stk::mesh::put_field(*procRankField,
  //                     metaData->element_rank(),
  //                     metaData->universal_part()
  //                     );
  //stk::io::set_field_role(*procRankField,
  //                        Ioss::Field::MESH);

  // real part
  Teuchos::RCP<ScalarFieldType> psir_field =
    Teuchos::rcpFromRef( metaData->declare_field<ScalarFieldType>( "psi_R" ) );
  stk::mesh::put_field(*psir_field,
                       metaData->node_rank(),
                       metaData->universal_part());
  stk::io::set_field_role(*psir_field,
                          Ioss::Field::TRANSIENT);

  // imaginary part
  Teuchos::RCP<ScalarFieldType> psii_field =
    Teuchos::rcpFromRef( metaData->declare_field<ScalarFieldType>( "psi_Z" ) );
  stk::mesh::put_field(*psii_field,
                       metaData->node_rank(),
                       metaData->universal_part());
  stk::io::set_field_role(*psii_field,
                          Ioss::Field::TRANSIENT);

  // Magnetic vector potential.
  // Think about whether to declare the magnetic vector potential as
  // Ioss::Field::TRANSIENT or Ioss::Field::ATTRIBUTE. "TRANSIENT" writes the
  // data out for all time/continuation steps (although) is it actually the
  // same throughout. "ATTRIBUTE" stores the vector field only once and hence
  // saves a lot of disk space, but not all readers recoginize the ATTRIBUTE
  // field (e.g., ParaView).
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
    Teuchos::rcpFromRef( metaData->declare_field<VectorFieldType>( "A" ) );
  stk::mesh::put_field(*mvpField,
                       metaData->node_rank(),
                       metaData->universal_part(),
                       numDim);
  stk::io::set_field_role(*mvpField,
                          Ioss::Field::ATTRIBUTE);

  // Sometimes, it may be stored in the file with three components (A_X, A_Y,
  // A_Z), sometimes, if the domain is two-dimensional, with two components
  // (typically A_R, A_Z then, for some funny reason -- cylindrical
  // coordinates?).
  //
  // To be on the safe side, declare the vector field A, and the scalar fields
  // A_R, A_Z here. Then further below apply some logic to make sense of the
  // findings.
  // What we'd really need here though is a "read everything that's in file"
  // kind of routine.
  Teuchos::RCP<ScalarFieldType> mvpFieldR =
    Teuchos::rcpFromRef( metaData->declare_field< ScalarFieldType >( "A_R" ) );
  stk::mesh::put_field(*mvpFieldR,
                       metaData->node_rank(),
                       metaData->universal_part());
  stk::io::set_field_role(*mvpFieldR,
                          Ioss::Field::ATTRIBUTE );

  Teuchos::RCP<ScalarFieldType> mvpFieldZ =
    Teuchos::rcpFromRef( metaData->declare_field< ScalarFieldType >( "A_Z" ) );
  stk::mesh::put_field(*mvpFieldZ,
                       metaData->node_rank(),
                       metaData->universal_part());
  stk::io::set_field_role(*mvpFieldZ,
                          Ioss::Field::ATTRIBUTE);

  // Thickness field. Same as above.
  Teuchos::RCP<ScalarFieldType> thicknessField =
    Teuchos::rcpFromRef(metaData->declare_field<ScalarFieldType>("thickness"));
  stk::mesh::put_field(*thicknessField,
                       metaData->node_rank(),
                       metaData->universal_part());
  stk::io::set_field_role(*thicknessField, Ioss::Field::ATTRIBUTE);

  // Potential field. Same as above.
  Teuchos::RCP<ScalarFieldType> potentialField =
    Teuchos::rcpFromRef(metaData->declare_field<ScalarFieldType>("V"));
  stk::mesh::put_field(*potentialField,
                       metaData->node_rank(),
                       metaData->universal_part());
  stk::io::set_field_role(*potentialField, Ioss::Field::ATTRIBUTE);
  // ---------------------------------------------------------------------------

  // define_input_fields() doesn't like the ATTRIBUTE fields; disable.
  // What was it good for anyways?
//  stk::io::define_input_fields( *meshData,
//                                *metaData
//                              );

//  stk::io::put_io_part_attribute( metaData->universal_part() );

  metaData->commit();

#ifdef HAVE_MPI
  if (fileIsSerial && comm.NumProc()>1)
  {
    bulkData->modification_begin();
    Ioss::Region *region = meshData->m_input_region;
    // Populate on the reader process.
    if (comm.MyPID() == 0)
      my_populate_bulk_data_(*bulkData, *region, *metaData);
    // Note:
    // Restart from a single Exodus file not currently supported.
    bulkData->modification_end();
  }
  else
  {
#endif
    stk::io::populate_bulk_data(*bulkData, *meshData);

    // Remember: Indices in STK are 1-based. :/
    stk::io::process_input_request(*meshData, *bulkData, index_+1);
    bulkData->modification_end();
#ifdef HAVE_MPI
  }

  // Rebalance.
  stk::mesh::Selector selector(metaData->universal_part());
  stk::mesh::Selector owned_selector(metaData->locally_owned_part());
  double imbalance = stk::rebalance::check_balance(*bulkData,
                                                   NULL,
                                                   metaData->node_rank(),
                                                   &selector);
  *out_ << "The imbalance is " << imbalance << ". ";
  if (imbalance < 1.5)
    *out_ << "That's acceptable." << std::endl;
  else
  {
    *out_ << "Rebalance!" << std::endl;
    // Zoltan graph-based reblancing.
    // http://trilinos.sandia.gov/packages/docs/dev/packages/stk/doc/html/group__stk__rebalance__unit__test__module.html
    Teuchos::ParameterList lb_method;
    lb_method.set("LOAD BALANCING METHOD", "4");
    Teuchos::ParameterList graph;
    graph.sublist(stk::rebalance::Zoltan::default_parameters_name()) = lb_method;
    stk::rebalance::Zoltan zoltan_partition(mcomm, numDim, graph);

    stk::rebalance::rebalance(*bulkData,
                              owned_selector,
                              &*coordinatesField,
                              NULL,
                              zoltan_partition);

    imbalance = stk::rebalance::check_balance(*bulkData,
                                              NULL,
                                              metaData->node_rank(),
                                              &selector);

    *out_ << "After rebalancing, the imbalance is " << imbalance << "." << std::endl;
  }
#endif

  // create the mesh with these specifications
  Teuchos::RCP<Nosh::StkMesh> mesh =
    Teuchos::rcp(new Nosh::StkMesh(comm, metaData, bulkData, coordinatesField));

  data.setName( "data" );

  // add the mesh to the data list
  data.set( "mesh", mesh );

  // create the state
  data.set("psi", this->complexfield2vector_(mesh, psir_field, psii_field));

  Teuchos::RCP<const Epetra_MultiVector> mvp;
  mvp = this->createMvp_(mesh, mvpField);
  // Check if it's 0.
  // If the field appears to be zeroed-out, it's probably not there.
  // Try A_R, A_Z.
  double r[3];
  TEUCHOS_ASSERT_EQUALITY(0, mvp->NormInf( r ));
  double tol = 1.0e-15;
  if (r[0]<tol && r[1]<tol && r[2]<tol)
    mvp = this->createMvpRZ_(mesh, mvpFieldR, mvpFieldZ);
  data.set( "A", mvp );

  // Check of the thickness data is of any value. If not: ditch it.
  Teuchos::RCP<Epetra_Vector> thickness = this->scalarfield2vector_(mesh, thicknessField);
  double norminf;
  TEUCHOS_ASSERT_EQUALITY(0, thickness->NormInf( &norminf ));
  if ( norminf < 1.0e-15 ) // assume that thickness wasn't present, fill with default value
    TEUCHOS_ASSERT_EQUALITY(0, thickness->PutScalar( 1.0 ));
  data.set( "thickness", thickness );
  // These are vain attempts to find out whether thicknessField is actually empty.
//     const stk::mesh::FieldBase::RestrictionVector & restrictions = thicknessField->restrictions();
//     TEUCHOS_ASSERT( !restrictions.empty() );
//     *out << "max_size " << thicknessField->max_size(metaData->node_rank()) << std::endl;

  // Check of the data is of any value. If not: ditch it.
  Teuchos::RCP<Epetra_Vector> potential = this->scalarfield2vector_(mesh, potentialField);
  TEUCHOS_ASSERT_EQUALITY(0, potential->NormInf( &norminf ));
  if ( norminf < 1.0e-15 ) // assume that potential wasn't present, fill with default value
    TEUCHOS_ASSERT_EQUALITY(0, potential->PutScalar( 0.0 ));
  data.set( "V", potential );

  return;
}
// =============================================================================
void
StkMeshReader::
my_populate_bulk_data_(stk::mesh::BulkData &bulk,
                       Ioss::Region &region,
                       stk::mesh::fem::FEMMetaData &fem_meta)
{
  // From Albany, Albany_IossSTKMeshStruct.cpp.
  // This function duplicates the function stk::io::populate_bulk_data, with the exception of an
  // internal modification_begin() and modification_end(). When reading bulk data on a single processor (single
  // exodus file), all PEs must enter the modification_begin() / modification_end() block, but only one reads the
  // bulk data.
  // TODO pull the modification statements from populate_bulk_data and retrofit.
  { // element blocks

    const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
    for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
  it != elem_blocks.end(); ++it) {
      Ioss::ElementBlock *entity = *it;

      if (stk::io::include_entity(entity)) {
  const std::string &name = entity->name();
  stk::mesh::Part* const part = fem_meta.get_part(name);
  assert(part != NULL);

  const CellTopologyData* cell_topo = stk::io::get_cell_topology(*part);
  if (cell_topo == NULL) {
    std::ostringstream msg ;
    msg << " INTERNAL_ERROR: Part " << part->name() << " returned NULL from get_cell_topology()";
    throw std::runtime_error( msg.str() );
  }

  std::vector<int> elem_ids ;
  std::vector<int> connectivity ;

  entity->get_field_data("ids", elem_ids);
  entity->get_field_data("connectivity", connectivity);

  size_t element_count = elem_ids.size();
  int nodes_per_elem = cell_topo->node_count ;

  std::vector<stk::mesh::EntityId> id_vec(nodes_per_elem);
  std::vector<stk::mesh::Entity*> elements(element_count);

  for(size_t i=0; i<element_count; ++i) {
    int *conn = &connectivity[i*nodes_per_elem];
    std::copy(&conn[0], &conn[0+nodes_per_elem], id_vec.begin());
    elements[i] = &stk::mesh::fem::declare_element(bulk, *part, elem_ids[i], &id_vec[0]);
  }

  // Add all element attributes as fields.
  // If the only attribute is 'attribute', then add it; otherwise the other attributes are the
  // named components of the 'attribute' field, so add them instead.
  Ioss::NameList names;
  entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
  for(Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
    if(*I == "attribute" && names.size() > 1)
      continue;
    stk::mesh::FieldBase *field = fem_meta.get_field<stk::mesh::FieldBase> (*I);
    if (field)
      stk::io::field_data_from_ioss(field, elements, entity, *I);
  }
      }
    }
  }

  { // nodeblocks

    const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
    assert(node_blocks.size() == 1);

    Ioss::NodeBlock *nb = node_blocks[0];

    std::vector<stk::mesh::Entity*> nodes;
    stk::io::get_entity_list(nb, fem_meta.node_rank(), bulk, nodes);

    stk::mesh::Field<double,stk::mesh::Cartesian> *coord_field =
      fem_meta.get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");

    stk::io::field_data_from_ioss(coord_field, nodes, nb, "mesh_model_coordinates");

  }

  { // nodesets

    const Ioss::NodeSetContainer& node_sets = region.get_nodesets();

    for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
  it != node_sets.end(); ++it) {
      Ioss::NodeSet *entity = *it;

      if (stk::io::include_entity(entity)) {
  const std::string & name = entity->name();
  stk::mesh::Part* const part = fem_meta.get_part(name);
  assert(part != NULL);
  stk::mesh::PartVector add_parts( 1 , part );

  std::vector<int> node_ids ;
  int node_count = entity->get_field_data("ids", node_ids);

  std::vector<stk::mesh::Entity*> nodes(node_count);
  stk::mesh::EntityRank n_rank = fem_meta.node_rank();
  for(int i=0; i<node_count; ++i) {
    nodes[i] = bulk.get_entity(n_rank, node_ids[i] );
    if (nodes[i] != NULL)
      bulk.declare_entity(n_rank, node_ids[i], add_parts );
  }

  stk::mesh::Field<double> *df_field =
    fem_meta.get_field<stk::mesh::Field<double> >("distribution_factors");

  if (df_field != NULL) {
    stk::io::field_data_from_ioss(df_field, nodes, entity, "distribution_factors");
  }
      }
    }
  }

  { // sidesets

    const Ioss::SideSetContainer& side_sets = region.get_sidesets();

    for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
  it != side_sets.end(); ++it) {
      Ioss::SideSet *entity = *it;

      if (stk::io::include_entity(entity)) {
//  process_surface_entity(entity, bulk);
  {
    assert(entity->type() == Ioss::SIDESET);

    size_t block_count = entity->block_count();
    for (size_t i=0; i < block_count; i++) {
      Ioss::SideBlock *block = entity->get_block(i);
      if (stk::io::include_entity(block)) {
  std::vector<int> side_ids ;
  std::vector<int> elem_side ;

  stk::mesh::Part * const sb_part = fem_meta.get_part(block->name());
  stk::mesh::EntityRank elem_rank = fem_meta.element_rank();

  block->get_field_data("ids", side_ids);
  block->get_field_data("element_side", elem_side);

  assert(side_ids.size() * 2 == elem_side.size());
  stk::mesh::PartVector add_parts( 1 , sb_part );

  size_t side_count = side_ids.size();
  std::vector<stk::mesh::Entity*> sides(side_count);
  for(size_t is=0; is<side_count; ++is) {
    stk::mesh::Entity* const elem = bulk.get_entity(elem_rank, elem_side[is*2]);

    // If NULL, then the element was probably assigned to an
    // element block that appears in the database, but was
    // subsetted out of the analysis mesh. Only process if
    // non-null.
    if (elem != NULL) {
      // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
      int side_ordinal = elem_side[is*2+1] - 1;

      stk::mesh::Entity* side_ptr = NULL;
      side_ptr = &stk::mesh::fem::declare_element_side(bulk, side_ids[is], *elem, side_ordinal);
      stk::mesh::Entity& side = *side_ptr;

      bulk.change_entity_parts( side, add_parts );
      sides[is] = &side;
    } else {
      sides[is] = NULL;
    }
  }

  const stk::mesh::Field<double, stk::mesh::ElementNode> *df_field =
    stk::io::get_distribution_factor_field(*sb_part);
  if (df_field != NULL) {
    stk::io::field_data_from_ioss(df_field, sides, block, "distribution_factors");
  }
      }
    }
  }
      }
    }
  }

}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
StkMeshReader::
complexfield2vector_( const Teuchos::RCP<const Nosh::StkMesh> &mesh,
                      const Teuchos::RCP<ScalarFieldType> &realField,
                      const Teuchos::RCP<ScalarFieldType> &imagField
                      ) const
{
  // Psi needs to have unique node IDs to be able to compute Norm2().
  // This is required in Belos.
  const std::vector<stk::mesh::Entity*> &ownedNodes = mesh->getOwnedNodes();

  // Create vector with this respective map.
  Teuchos::RCP<Epetra_Vector> vector =
    Teuchos::rcp( new Epetra_Vector( *mesh->getComplexNonOverlapMap() ) );

  // Fill the vector with data from the file.
  for ( unsigned int k=0; k<ownedNodes.size(); k++ )
  {
    // real part
    double* realVal = stk::mesh::field_data(*realField, *ownedNodes[k]);
    (*vector)[2*k] = realVal[0];

    // imaginary part
    double* imagVal = stk::mesh::field_data( *imagField, *ownedNodes[k] );
    (*vector)[2*k+1] = imagVal[0];
  }

#ifdef _DEBUG_
  double r;
  TEUCHOS_ASSERT_EQUALITY(0, vector->NormInf( &r ));
  TEUCHOS_TEST_FOR_EXCEPT_MSG( r!=r || r>1.0e100,
                       "The input data seems flawed. Abort." );
#endif

  return vector;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
StkMeshReader::
scalarfield2vector_( const Teuchos::RCP<const Nosh::StkMesh> &mesh,
                     const Teuchos::RCP<ScalarFieldType> &field
                    ) const
{
  // Get overlap nodes.
  const std::vector<stk::mesh::Entity*> &overlapNodes = mesh->getOverlapNodes();

  // Create vector with this respective map.
  Teuchos::RCP<Epetra_Vector> vector =
    Teuchos::rcp( new Epetra_Vector( *mesh->getNodesOverlapMap() ) );

#ifdef _DEBUG_
  TEUCHOS_ASSERT( !field.is_null() );
#endif

  // Fill the vector with data from the file.
  for ( unsigned int k=0; k<overlapNodes.size(); k++ )
  {
    double* fieldVal = stk::mesh::field_data( *field,
                                              *overlapNodes[k] );
    // Check if the field is actually there.
    if (fieldVal == NULL)
    {
      *out_ << "WARNING: Thickness value for node " << k << " not found.\n"
      <<
      "Probably there is no field given with the state. Using default."
      << std::endl;
      return Teuchos::null;
    }
    (*vector)[k] = fieldVal[0];
  }

#ifdef _DEBUG_
  double r;
  // Use NormInf as it's robust against overlapping maps.
  TEUCHOS_ASSERT_EQUALITY(0, vector->NormInf( &r ));
  TEUCHOS_TEST_FOR_EXCEPT_MSG( r!=r || r>1.0e100,
                       "The input data seems flawed. Abort." );
#endif

  return vector;
}
// =============================================================================
Teuchos::RCP<Epetra_MultiVector>
StkMeshReader::
createMvp_( const Teuchos::RCP<const Nosh::StkMesh> &mesh,
            const Teuchos::RCP<const VectorFieldType> &mvpField
            ) const
{
  // Get overlap nodes.
  const std::vector<stk::mesh::Entity*> &overlapNodes = mesh->getOverlapNodes();

  // Create vector with this respective map.
  int numComponents = 3;
  Teuchos::RCP<Epetra_MultiVector> mvp =
    Teuchos::rcp( new Epetra_MultiVector( *mesh->getNodesOverlapMap(),
                                          numComponents ) );

#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mvpField.is_null() );
#endif

  // Fill the vector with data from the file.
  for ( unsigned int k=0; k<overlapNodes.size(); k++ )
  {
    const double * const mvpVal =
      stk::mesh::field_data( *mvpField, *overlapNodes[k] );
#ifdef _DEBUG_
    // Check if the field is actually there.
    TEUCHOS_TEST_FOR_EXCEPT_MSG( mvpVal == NULL,
      "MVPX value for node " << k << " not found.\n" <<
      "Probably there is no MVP field given with the state."
      );
#endif
    (*mvp)[0][k] = mvpVal[0];
    (*mvp)[1][k] = mvpVal[1];
    (*mvp)[2][k] = mvpVal[2];
  }

#ifdef _DEBUG_
  // Check for NaNs and uninitialized data.
  double r[3];
  // Use NormInf as it's robust against overlapping maps.
  TEUCHOS_ASSERT_EQUALITY(0, mvp->NormInf( r ));
  TEUCHOS_TEST_FOR_EXCEPT_MSG( r[0]!=r[0] || r[0]>1.0e100
                       || r[1]!=r[1] || r[1]>1.0e100
                       || r[2]!=r[2] || r[2]>1.0e100,
                       "The input data seems flawed. Abort." );
#endif

  return mvp;
}
// =============================================================================
Teuchos::RCP<Epetra_MultiVector>
StkMeshReader::
createMvpRZ_( const Teuchos::RCP<const Nosh::StkMesh> &mesh,
              const Teuchos::RCP<const ScalarFieldType> &mvpFieldR,
              const Teuchos::RCP<const ScalarFieldType> &mvpFieldZ
            ) const
{
  // Get overlap nodes.
  const std::vector<stk::mesh::Entity*> &overlapNodes = mesh->getOverlapNodes();

  // Create vector with this respective map.
  int numComponents = 3;
  Teuchos::RCP<Epetra_MultiVector> mvp =
    Teuchos::rcp( new Epetra_MultiVector( *mesh->getNodesOverlapMap(),
                                          numComponents ) );

#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mvpFieldR.is_null() );
  TEUCHOS_ASSERT( !mvpFieldZ.is_null() );
#endif
  // Fill the vector with data from the file.
  for ( unsigned int k=0; k<overlapNodes.size(); k++ )
  {
    double *mvpValR = stk::mesh::field_data( *mvpFieldR, *overlapNodes[k] );
    double *mvpValZ = stk::mesh::field_data( *mvpFieldZ, *overlapNodes[k] );
#ifdef _DEBUG_
    // Check if the field is actually there.
    TEUCHOS_TEST_FOR_EXCEPT_MSG( mvpValR == NULL,
      "MVPR value for node " << k << " not found.\n" <<
      "Probably there is no MVP field given with the state."
      );
    TEUCHOS_TEST_FOR_EXCEPT_MSG( mvpValZ == NULL,
      "MVPZ value for node " << k << " not found.\n" <<
      "Probably there is no MVP field given with the state."
      );
#endif
    (*mvp)[0][k] = *mvpValR;
    (*mvp)[1][k] = *mvpValR;
    // TODO Remove explicit extension by 0.
    (*mvp)[2][k] = 0.0;
  }

#ifdef _DEBUG_
  // Check for NaNs and uninitialized data.
  double r[2];
  // Use NormInf as it's robust against overlapping maps.
  TEUCHOS_ASSERT_EQUALITY(0, mvp->NormInf( r ));
  TEUCHOS_TEST_FOR_EXCEPT_MSG( r[0]!=r[0] || r[0]>1.0e100
                       || r[1]!=r[1] || r[1]>1.0e100,
                       "The input data seems flawed. Abort." );
#endif

  return mvp;
}
// =============================================================================
} // namespace Nosh
// =============================================================================
//
// Helper functions
//
void
Nosh::
StkMeshRead( const Epetra_Comm &comm,
             const std::string &fileName,
             const int step,
             Teuchos::ParameterList &data
             )
{
  StkMeshReader reader( fileName, step);
  reader.read(comm, data);
  return;
}
// =============================================================================
