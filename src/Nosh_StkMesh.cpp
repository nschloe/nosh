// @HEADER
//
//    Mesh class with compatibility to stk_mesh.
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
#include "Nosh_StkMesh.hpp"

#include <Trilinos_version.h>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Export.h>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialSpdDenseSolver.hpp>

// #include <stk_mesh/fem/FEMMetaData.hpp>
// #include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
// #include <stk_mesh/base/Comm.hpp> // for comm_mesh_counts
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>
#include <stk_mesh/base/GetEntities.hpp>
// For parallel_reduce and sum:
#include <stk_mesh/base/FieldParallel.hpp>

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

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  #include <Teuchos_TimeMonitor.hpp>
#endif
// =============================================================================
namespace Nosh {
// =============================================================================
StkMesh::
StkMesh(const Epetra_Comm & comm,
        const std::string & fileName,
        const int index
        ) :
  numDim_( 3 ),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  computeEdgeCoefficientsTime_(Teuchos::TimeMonitor::getNewTimer(
                                  "Nosh: StkMesh::computeEdgeCoefficients")),
  writeTime_(Teuchos::TimeMonitor::getNewTimer(
                                  "Nosh: StkMesh::write")),
#endif
  comm_( comm ),
  metaData_( numDim_ ),
  meshData_(Teuchos::rcp(new stk::io::MeshData())),
  bulkData_(stk::mesh::fem::FEMMetaData::get_meta_data(metaData_),
#ifdef HAVE_MPI
    Teuchos::dyn_cast<const Epetra_MpiComm>(comm).Comm()
#else
    1
#endif
  ),
  ownedNodes_(),
  nodesMap_(),
  nodesOverlapMap_(),
  complexMap_(),
  complexOverlapMap_(),
  fvmEntitiesUpToDate_( false ),
  controlVolumes_(),
  controlVolumesUpToDate_( false ),
  edgeCoefficients_(),
  edgeCoefficientsUpToDate_( false ),
  edgeNodes_(),
  cellEdges_( Teuchos::null ),
  outputChannelIsOpen_(false),
  time_(0.0)
{
  this->read(comm, fileName, index);
  // Create all the maps.
  ownedNodes_ = buildOwnedNodes_();
  nodesMap_ = this->createEntitiesMap_(ownedNodes_);
  nodesOverlapMap_ = this->createEntitiesMap_(this->getOverlapNodes());
  complexMap_ = this->createComplexMap_(ownedNodes_);
  complexOverlapMap_ = this->createComplexMap_( this->getOverlapNodes() );
  controlVolumes_ = Teuchos::rcp( new Epetra_Vector( *nodesMap_ ) );
  // Create adjacent entities.
  this->createEdges_();
//  int nodesPerCell;
//  if (comm_.MyPID() == 0)
//  {
//    std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
//    nodesPerCell = cells[0]->relations(metaData_.node_rank()).size();
//  }
//  comm_.Broadcast(&nodesPerCell, 1, 0);
//  if (nodesPerCell >= 4)
//  {
//    stk::mesh::PartVector add_parts;
//    stk::mesh::create_adjacent_entities(bulkData_, add_parts);
//  }
}
// =============================================================================
StkMesh::
~StkMesh()
{
}
// =============================================================================
void
StkMesh::
read(const Epetra_Comm &comm,
     const std::string &fileName,
     const int index
     )
{
#ifdef HAVE_MPI
  const Epetra_MpiComm &mpicomm =
    Teuchos::dyn_cast<const Epetra_MpiComm>(comm);
  MPI_Comm mcomm = mpicomm.Comm();
#endif
  // Take two different fields with one component
  // instead of one field with two components. This works around
  // Ioss's inability to properly read psi_R, psi_Z as a complex variable.
  // (It can handle data_X, data_Y, data_Z though.)
  //const unsigned int neq = 1;

  // ---------------------------------------------------------------------------
  // initialize database communication
  Ioss::Init::Initializer io;

  // If the file is serial, read it with process 0 and embed it
  // in the multiproc context. Load balancing is done later anyways.
  MPI_Comm readerComm = bulkData_.parallel();
#ifdef HAVE_MPI
  bool fileIsSerial = fileName.substr(fileName.find_last_of(".") + 1) == "e";
  if (fileIsSerial && comm.NumProc() > 1)
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

  std::string meshType = "exodusii";
  // By default, the Exodus field separator is '_'
  // such that only fields 'A_*' would be recognized as such.
  // However, VTK is inconsistent in its putting underscores
  // in field names, such that vectors with three components
  // are typically stored as 'AX', 'AY', 'AZ'.
  // This can be worked around where the Exodus file is created
  // or here.
  // The problem about messing around with the input database here
  // is that the output routines will always "correctly" append
  // an underscore. This then creates a difference between VTK-
  // generated and Trilinos-generated files.
  // Setting removing the underscore from both input and output
  // files is possible in a figure version of Trilinos, c.f.
  // Greg's mail (2012-09-20):
  // '''
  // Basically, if you do the following prior to creating ,
  // it should give you a consistent field separator:
  //
  // meshData_->m_property_manager.add(Ioss::Property("FIELD_SUFFIX_SEPARATOR", ""));
  // '''
  //
  // Anyways. Here's some code that removes the underscore from
  // the input database. Keep it commented out for now though.
  // <WORKAROUND CODE START>
  //Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(meshType,
  //                                                fileName,
  //                                                Ioss::READ_MODEL,
  //                                                readerComm);
  //TEUCHOS_TEST_FOR_EXCEPT_MSG( dbi == NULL || !dbi->ok(),
  //                     "ERROR: Could not open database '" << fileName
  //                     << "' of type '" << meshType << "'.");
  //// set the vector field label separator
  //dbi->set_field_separator( 0 );
  //// create region to feed into meshData
  //Ioss::Region *in_region = new Ioss::Region( dbi, "input_model" );
  //meshData_->m_input_region = in_region;
  // <WORKAROUND CODE END>

  // This checks the existence of the file, checks to see if we can open it,
  // builds a handle to the region and puts it in mesh_data (in_region),
  // and reads the metaData into metaData.
  stk::io::create_input_mesh(meshType,
                             fileName,
                             readerComm,
                             metaData_,
                             *meshData_);
  // ---------------------------------------------------------------------------
  // As of now (2012-03-21) there is no way to determine which fields are
  // actually present in the file. The only thing we can do is to declare them,
  // and check if they're full of zeros in the end.
  VectorFieldType &coordinatesField =
    metaData_.declare_field<VectorFieldType>("coordinates");
  stk::mesh::put_field(coordinatesField,
                       metaData_.node_rank(),
                       metaData_.universal_part(),
                       numDim_);
  stk::io::set_field_role(coordinatesField,
                          Ioss::Field::MESH);

  IntScalarFieldType &procRankField =
    metaData_.declare_field<IntScalarFieldType>("proc_rank");
  stk::mesh::put_field(procRankField,
                       metaData_.element_rank(),
                       metaData_.universal_part()
                       );
  stk::io::set_field_role(procRankField,
                          Ioss::Field::MESH);

  // real part
  ScalarFieldType &psir_field =
    metaData_.declare_field<ScalarFieldType>("psi_R");
  stk::mesh::put_field(psir_field,
                       metaData_.node_rank(),
                       metaData_.universal_part());
  stk::io::set_field_role(psir_field,
                          Ioss::Field::TRANSIENT);

  // imaginary part
  ScalarFieldType &psii_field =
    metaData_.declare_field<ScalarFieldType>("psi_Z");
  stk::mesh::put_field(psii_field,
                       metaData_.node_rank(),
                       metaData_.universal_part());
  stk::io::set_field_role(psii_field,
                          Ioss::Field::TRANSIENT);

  // Magnetic vector potential.
  // Unconditionally assume that the field is 3D (A_X, A_Y, A_Z) even
  // if the domain is two-dimensional. Eventually, we only need the
  // the projection of the field onto the edges of the mesh, this
  // this might be a bit overkill. Until we can explicitly associate
  // fields with edges, though, keep it this way.
  // Also, declare "A" as Ioss::Field::ATTRIBUTE. This makes sure that
  // the data is written out, but only once (hence "attribute") and not
  // once per step. (This is with trilinos-dev as of July 2012.)
  VectorFieldType &mvpField =
    metaData_.declare_field<VectorFieldType>("A");
  stk::mesh::put_field(mvpField,
                       metaData_.node_rank(),
                       metaData_.universal_part(),
                       numDim_);
  stk::io::set_field_role(mvpField,
                          Ioss::Field::ATTRIBUTE);

  // Note:
  // Thickness and V are actually scalar fields. However, they are both
  // declared VectorFields with 1 component here since otherwise
  // SEACAS throws
  //
  //   p=0: *** Caught standard std::exception of type 'std::runtime_error' :
  //
  //   Expr '!(restr->not_equal_stride(tmp))' eval'd to true, throwing.
  //   Error occured at: stk_mesh/stk_mesh/baseImpl/FieldBaseImpl.cpp:240
  //   Error: std::mesh::MetaData::declare_field_restriction FAILED for FieldBaseImpl<double,Cartesian3d>[ name = "thickness" , #states = 1 ] { entity_rank(0) part({UNIVERSAL}) : 1 } WITH INCOMPATIBLE REDECLARATION { entity_rank(0) part({UNIVERSAL}) : 19754528 }
  //
  // Thickness field. Same as above.
  VectorFieldType &thicknessField =
    metaData_.declare_field<VectorFieldType>("thickness");
  stk::mesh::put_field(thicknessField,
                       metaData_.node_rank(),
                       metaData_.universal_part(),
                       1);
  stk::io::set_field_role(thicknessField, Ioss::Field::ATTRIBUTE);

  // Potential field. Same as above.
  VectorFieldType &potentialField =
    metaData_.declare_field<VectorFieldType>("V");
  stk::mesh::put_field(potentialField,
                       metaData_.node_rank(),
                       metaData_.universal_part(),
                       1);
  stk::io::set_field_role(potentialField, Ioss::Field::ATTRIBUTE);
  // ---------------------------------------------------------------------------

  // define_input_fields() doesn't like the ATTRIBUTE fields; disable.
  // What was it good for anyways?
//  stk::io::define_input_fields( *meshData_,
//                                metaData
//                              );

//  stk::io::put_io_part_attribute( metaData_.universal_part() );

  // Finalize the setup.
  metaData_.commit();

#ifdef HAVE_MPI
  if (fileIsSerial && comm.NumProc() > 1)
  {
    bulkData_.modification_begin();
    Ioss::Region *region = meshData_->m_input_region;
    if (comm.MyPID() == 0)
    {
      stk::io::process_mesh_bulk_data(region, bulkData_);
      stk::io::input_mesh_fields(region, bulkData_, index+1);
    }
    bulkData_.modification_end();
  }
  else
  {
#endif
    stk::io::populate_bulk_data(bulkData_, *meshData_);

    // Remember: Indices in STK are 1-based. :/
    stk::io::process_input_request(*meshData_, bulkData_, index+1);

    // This should be propagated into stk::io
    //Ioss::Region *region = meshData_->m_input_region;
    //const Ioss::ElementBlockContainer& elem_blocks = region->get_element_blocks();
    //// Uncomment to print what fields are in the exodus file
    //Ioss::NameList exo_fld_names;
    //elem_blocks[0]->field_describe(&exo_fld_names);
    //for(std::size_t i = 0; i < exo_fld_names.size(); i++)
    //  std::cout << "Found field \"" << exo_fld_names[i] << "\" in exodus file" << std::endl;

    bulkData_.modification_end();
#ifdef HAVE_MPI
  }

  // Rebalance.
  stk::mesh::Selector selector(metaData_.universal_part());
  stk::mesh::Selector owned_selector(metaData_.locally_owned_part());
  double imbalance = stk::rebalance::check_balance(bulkData_,
                                                   NULL,
                                                   metaData_.node_rank(),
                                                   &selector);
  if (imbalance > 1.5)
  {
    //if (comm.MyPID() == 0)
    //  std::cout << "The imbalance is " << imbalance << ". Rebalance!" << std::endl;
    // Zoltan graph-based reblancing.
    // http://trilinos.sandia.gov/packages/docs/dev/packages/stk/doc/html/group__stk__rebalance__unit__test__module.html
    Teuchos::ParameterList lb_method;
    lb_method.set("LOAD BALANCING METHOD", "4");
    Teuchos::ParameterList graph;
    graph.sublist(stk::rebalance::Zoltan::default_parameters_name()) = lb_method;
    stk::rebalance::Zoltan zoltan_partition(mcomm, numDim_, graph);

    stk::rebalance::rebalance(bulkData_,
                              owned_selector,
                              &coordinatesField,
                              NULL,
                              zoltan_partition);

    //imbalance = stk::rebalance::check_balance(bulkData_,
    //                                          NULL,
    //                                          metaData_.node_rank(),
    //                                          &selector);
    //if (comm.MyPID() == 0)
    //  std::cout << "After rebalancing, the imbalance is " << imbalance << "." << std::endl;
  }
#endif

  // Extract time value.
  time_ = meshData_->m_input_region->get_state_time(index+1);

  return;
}
// =============================================================================
double
StkMesh::
getTime() const
{
  return time_;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
StkMesh::
complexfield2vector_(const ScalarFieldType &realField,
                     const ScalarFieldType &imagField
                     ) const
{
  // Psi needs to have unique node IDs to be able to compute Norm2().
  // This is required in Belos.
  const std::vector<stk::mesh::Entity*> &ownedNodes = this->getOwnedNodes();

  // Create vector with this respective map.
  Teuchos::RCP<Epetra_Vector> vector =
    Teuchos::rcp(new Epetra_Vector(*this->getComplexNonOverlapMap()));

  // Fill the vector with data from the file.
  for (unsigned int k=0; k<ownedNodes.size(); k++)
  {
    // real part
    double* realVal = stk::mesh::field_data(realField, *ownedNodes[k]);
    (*vector)[2*k] = realVal[0];

    // imaginary part
    double* imagVal = stk::mesh::field_data(imagField, *ownedNodes[k]);
    (*vector)[2*k+1] = imagVal[0];
  }

#ifndef NDEBUG
  double r;
  TEUCHOS_ASSERT_EQUALITY(0, vector->NormInf( &r ));
  TEUCHOS_TEST_FOR_EXCEPT_MSG( r!=r || r>1.0e100,
                       "The input data seems flawed. Abort." );
#endif

  return vector;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
StkMesh::
field2vector_(const ScalarFieldType &field) const
{
  // Get overlap nodes.
  const std::vector<stk::mesh::Entity*> &overlapNodes =
    this->getOverlapNodes();

  // Create vector with this respective map.
  Teuchos::RCP<Epetra_Vector> vector =
    Teuchos::rcp(new Epetra_Vector(*this->getNodesOverlapMap()));

  // Fill the vector with data from the file.
  for ( unsigned int k=0; k<overlapNodes.size(); k++ )
  {
    double* vals = stk::mesh::field_data(field,
                                         *overlapNodes[k] );
    // Check if the field is actually there.
#ifndef NDEBUG
    TEUCHOS_ASSERT(vals != NULL);
    //*out << "WARNING: Value for node " << k << " not found.\n" <<
    //  "Probably there is no field given with the state. Using default."
    //  << std::endl;
#endif
    (*vector)[k] = vals[0];
  }

#ifndef NDEBUG
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
StkMesh::
field2vector_(const VectorFieldType &field,
              const int numComponents
              ) const
{
  // Get overlap nodes.
  const std::vector<stk::mesh::Entity*> &overlapNodes =
    this->getOverlapNodes();

  // Create vector with this respective map.
  Teuchos::RCP<Epetra_MultiVector> vector =
    Teuchos::rcp(new Epetra_MultiVector(*(this->getNodesOverlapMap()),
                                         numComponents));

  // Fill the vector with data from the file.
  for ( unsigned int k=0; k<overlapNodes.size(); k++ )
  {
    const double * const vals =
      stk::mesh::field_data(field, *overlapNodes[k]);
#ifndef NDEBUG
    // Check if the field is actually there.
    TEUCHOS_TEST_FOR_EXCEPT_MSG( vals == NULL,
      "Field value for node " << k << " not found.\n" <<
      "Probably there is no field given with the state."
      );
#endif
    // Copy over.
    // A multivector isn't actually a good data structure for this.
    // What would be needed is a vector where each entry has k
    // components. This way, the data could stick together.
    for (int i=0; i<numComponents; i++)
      (*vector)[i][k] = vals[i];
  }

#ifndef NDEBUG
  // Check for NaNs and uninitialized data.
  std::vector<double> r(numComponents);
  // Use NormInf as it's robust against overlapping maps.
  TEUCHOS_ASSERT_EQUALITY(0, vector->NormInf(&r[0]));
  bool makesSense = true;
  for (int i=0; i<numComponents; i++)
  {
    if (r[i]!=r[i] || r[i]>1.0e100)
    {
      makesSense = false;
      break;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!makesSense,
                              "The input data seems flawed. Abort.");
#endif

  return vector;
}
// =============================================================================
void
StkMesh::
openOutputChannel(const string &outputDir,
                  const string &fileBaseName
                  )
{
  // prepare the data for output
#ifdef HAVE_MPI
  const Epetra_MpiComm &mpicomm =
    Teuchos::dyn_cast<const Epetra_MpiComm>(comm_);
  MPI_Comm mcomm = mpicomm.Comm();
#else
  const int mcomm = 1;
#endif
  const std::string extension = ".e";

  // Make sure the outputDir ends in "/".
  // Dir and filename are not concatenated properly in stk::mesh,
  std::stringstream outputFile;
  outputFile << outputDir << "/" << fileBaseName << extension;

  stk::io::create_output_mesh(outputFile.str(),
                              mcomm,
                              bulkData_,
                              *meshData_
                              );

  stk::io::define_output_fields(*meshData_,
                                metaData_
                                );

  outputChannelIsOpen_ = true;
  return;
}
// =============================================================================
void
StkMesh::
write(const Epetra_Vector & psi,
      const double time
      ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm( *writeTime_ );
#endif

  // Merge the state into the mesh.
//     mesh_->getBulkData()->modification_begin();
  this->mergeComplexVector_(psi, "psi");
//     mesh_->getBulkData()->modification_end();

  TEUCHOS_ASSERT(outputChannelIsOpen_);

  // Write it out to the file that's been specified in mesh_.
  const int out_step =
    stk::io::process_output_request(*meshData_, bulkData_, time);

  return;
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
StkMesh::
createVector(const std::string & fieldName) const
{
  const ScalarFieldType * const field =
    metaData_.get_field<ScalarFieldType>(fieldName);
#ifndef NDEBUG
  TEUCHOS_ASSERT(field != NULL);
#endif

  return this->field2vector_(*field);
}
// =============================================================================
Teuchos::RCP<Epetra_MultiVector>
StkMesh::
createMultiVector(const std::string & fieldName) const
{
  const VectorFieldType * const field =
    metaData_.get_field<VectorFieldType>(fieldName);
#ifndef NDEBUG
  TEUCHOS_ASSERT(field != NULL);
#endif

  return this->field2vector_(*field, 3);
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
StkMesh::
createComplexVector(const std::string & fieldName) const
{
  const ScalarFieldType * const r_field =
    metaData_.get_field<ScalarFieldType>(fieldName + "_R");
  const ScalarFieldType * const i_field =
    metaData_.get_field<ScalarFieldType>(fieldName + "_Z");
#ifndef NDEBUG
  TEUCHOS_ASSERT(r_field != NULL);
  TEUCHOS_ASSERT(i_field != NULL);
#endif

  return this->complexfield2vector_(*r_field, *i_field);
}
// =============================================================================
void
StkMesh::
mergeComplexVector_(const Epetra_Vector & psi,
                    const std::string & fieldName
                    ) const
{
  ScalarFieldType * psir_field =
    metaData_.get_field<ScalarFieldType>(fieldName + "_R");
  ScalarFieldType * psii_field =
    metaData_.get_field<ScalarFieldType>(fieldName + "_Z");
#ifndef NDEBUG
  TEUCHOS_ASSERT( psir_field != NULL );
  TEUCHOS_ASSERT( psii_field != NULL );
#endif

  // Zero out all nodal values, including the overlaps.
  const std::vector<stk::mesh::Entity*> &overlapNodes = this->getOverlapNodes();
  for (unsigned int k=0; k < overlapNodes.size(); k++)
  {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data( *psir_field, *overlapNodes[k] );
    *localPsiR = 0.0;
    double* localPsiI = stk::mesh::field_data( *psii_field, *overlapNodes[k] );
    *localPsiI = 0.0;
  }

  // Set owned nodes.
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY((unsigned int)psi.MyLength(), 2*ownedNodes_.size());
#endif
  for (unsigned int k=0; k < ownedNodes_.size(); k++)
  {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data( *psir_field, *ownedNodes_[k] );
    *localPsiR = psi[2*k];
    double* localPsiI = stk::mesh::field_data( *psii_field, *ownedNodes_[k] );
    *localPsiI = psi[2*k+1];
  }

  // This communication updates the field values on un-owned nodes
  // it is correct because the zeroSolutionField above zeros them all
  // and the getSolutionField only sets the owned nodes.
  stk::mesh::parallel_reduce(bulkData_,
                             stk::mesh::sum(*psir_field));
  stk::mesh::parallel_reduce(bulkData_,
                             stk::mesh::sum(*psii_field));

  return;
}
// =============================================================================
unsigned int
StkMesh::
getNumNodes() const
{
  return nodesMap_->NumGlobalElements();
}
// =============================================================================
Teuchos::RCP<const Epetra_Vector>
StkMesh::
getControlVolumes() const
{
  if ( !controlVolumesUpToDate_ )
    this->computeControlVolumes_();
  return controlVolumes_;
}
// =============================================================================
double
StkMesh::
getDomainVolume() const
{
  if ( !controlVolumesUpToDate_)
    this->computeControlVolumes_();
  // update the domain area value
  double volume;
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->Norm1( &volume ) );

  return volume;
}
// =============================================================================
const Epetra_Comm &
StkMesh::
getComm() const
{
  return comm_;
}
// =============================================================================
Teuchos::ArrayRCP<double>
StkMesh::
getEdgeCoefficients() const
{
  if ( !edgeCoefficientsUpToDate_ )
    this->computeEdgeCoefficients_();
  return edgeCoefficients_;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
StkMesh::
getOwnedCells() const
{
  // get owned elements
  stk::mesh::Selector select_owned_in_part =
      stk::mesh::Selector(metaData_.universal_part() )
    & stk::mesh::Selector(metaData_.locally_owned_part());
  std::vector<stk::mesh::Entity*> cells;
  stk::mesh::get_selected_entities(select_owned_in_part,
                                   bulkData_.buckets( metaData_.element_rank() ),
                                   cells
                                   );
  return cells;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
StkMesh::
getOverlapEdges() const
{
  // get overlap edges
  stk::mesh::Selector select_overlap_in_part =
       stk::mesh::Selector(metaData_.universal_part() )
    & (stk::mesh::Selector(metaData_.locally_owned_part())
      |stk::mesh::Selector(metaData_.globally_shared_part()));

  std::vector<stk::mesh::Entity*> edges;
  stk::mesh::get_selected_entities(select_overlap_in_part,
                                   bulkData_.buckets( metaData_.edge_rank() ),
                                   edges
                                   );
  return edges;
}
// =============================================================================
const Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> >
StkMesh::
getEdgeNodes() const
{
  return edgeNodes_;
}
// =============================================================================
double
StkMesh::
getScalarFieldNonconst(const stk::mesh::Entity * nodeEntity,
                       const std::string & fieldName
                       ) const
{
  const ScalarFieldType * const field =
    metaData_.get_field<ScalarFieldType>(fieldName);
#ifndef NDEBUG
  TEUCHOS_ASSERT(field != NULL);
#endif
  return *stk::mesh::field_data(*field, *nodeEntity);
}
// =============================================================================
const DoubleVector
StkMesh::
getVectorFieldNonconst(const stk::mesh::Entity * nodeEntity,
                       const std::string & fieldName,
                       const int numDims
                       ) const
{
  const VectorFieldType * const field =
    metaData_.get_field<VectorFieldType>(fieldName);
#ifndef NDEBUG
  TEUCHOS_ASSERT(field != NULL);
#endif
  // Make a Teuchos::Copy here as the access is nonconst.
  return DoubleVector(Teuchos::Copy,
                      stk::mesh::field_data(*field, *nodeEntity),
                      numDims);
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
StkMesh::
getNodesMap() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !nodesMap_.is_null() );
#endif
  return nodesMap_;
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
StkMesh::
getNodesOverlapMap() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !nodesOverlapMap_.is_null() );
#endif
  return nodesOverlapMap_;
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
StkMesh::
getComplexNonOverlapMap() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !complexMap_.is_null() );
#endif
  return complexMap_;
}
// =============================================================================
Teuchos::RCP<const Epetra_Map>
StkMesh::
getComplexOverlapMap() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !complexOverlapMap_.is_null() );
#endif
  return complexOverlapMap_;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
StkMesh::
buildOwnedNodes_() const
{
  stk::mesh::Selector select_owned_in_part =
        stk::mesh::Selector(metaData_.universal_part())
      & stk::mesh::Selector(metaData_.locally_owned_part());

  std::vector<stk::mesh::Entity*> on;
  stk::mesh::get_selected_entities(select_owned_in_part,
                                   bulkData_.buckets( metaData_.node_rank() ),
                                   on);
  return on;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
StkMesh::
getOwnedNodes() const
{
  return ownedNodes_;
}
// =============================================================================
std::vector<stk::mesh::Entity*>
StkMesh::
getOverlapNodes() const
{
  //  overlapnodes used for overlap map -- stored for changing coords
  stk::mesh::Selector select_overlap_in_part =
       stk::mesh::Selector(metaData_.universal_part())
    & (stk::mesh::Selector(metaData_.locally_owned_part())
      |stk::mesh::Selector(metaData_.globally_shared_part()));

  std::vector<stk::mesh::Entity*> overlapNodes;
  stk::mesh::get_selected_entities(select_overlap_in_part,
                                   bulkData_.buckets( metaData_.node_rank() ),
                                   overlapNodes
                                   );

  return overlapNodes;
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
StkMesh::
createEntitiesMap_( const std::vector<stk::mesh::Entity*> &entityList ) const
{
  const int numEntities = entityList.size();
  Teuchos::Array<int> gids(numEntities);
  for (int i=0; i < numEntities; i++)
    gids[i] = entityList[i]->identifier() - 1;

  return Teuchos::rcp(new Epetra_Map(-1, numEntities, gids.getRawPtr(), 0, comm_));
}
// =============================================================================
Teuchos::RCP<Epetra_Map>
StkMesh::
createComplexMap_( const std::vector<stk::mesh::Entity*> &nodeList ) const
{
  // Create a map for real/imaginary out of this.
  const int numDof = 2 * nodeList.size();
  Teuchos::Array<int> gids(numDof);
  for ( unsigned int k=0; k < nodeList.size(); k++ )
  {
    int globalNodeId = nodeList[k]->identifier() - 1;
    gids[2*k]   = 2*globalNodeId;
    gids[2*k+1] = 2*globalNodeId + 1;
  }
  return Teuchos::rcp(new Epetra_Map(-1, numDof, gids.getRawPtr(), 0, comm_));
}
// =============================================================================
unsigned int
StkMesh::
getNumEdgesPerCell( unsigned int cellDimension ) const
{
  // In n-simplices, all nodes are connected with all other nodesMap.
  // Hence, numEdges==sum_{i=1}^(numLocalNodes-1) i.
  unsigned int numEdgesPerCell = 0;
  for ( unsigned int i=1; i<cellDimension+1; i++ )
    numEdgesPerCell += i;

  return numEdgesPerCell;
}
// =============================================================================
void
StkMesh::
computeEdgeCoefficients_() const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm( *computeEdgeCoefficientsTime_ );
#endif

  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  Teuchos::Array<Teuchos::Tuple<stk::mesh::Entity*,2> >::size_type numEdges = edgeNodes_.size();

  if ( edgeCoefficients_.size() != numEdges )
    edgeCoefficients_.resize( numEdges );

  // Calculate the contributions edge by edge.
  for (unsigned int k=0; k < numCells; k++)
  {
    // Get edge coordinates.
    unsigned int numLocalEdges = cellEdges_[k].size();
    Teuchos::ArrayRCP<DoubleVector> localEdgeCoords( numLocalEdges );
    for ( unsigned int i=0; i<numLocalEdges; i++)
    {
      localEdgeCoords[i] = this->getVectorFieldNonconst(edgeNodes_[cellEdges_[k][i]][1],
                                                        "coordinates",
                                                        3);
      localEdgeCoords[i] -= this->getVectorFieldNonconst(edgeNodes_[cellEdges_[k][i]][0],
                                                         "coordinates",
                                                         3);
    }

    DoubleVector edgeCoeffs = getEdgeCoefficientsNumerically_( localEdgeCoords );

    // Fill the edge coefficients into the vector.
    for ( unsigned int i=0; i<numLocalEdges; i++ )
      edgeCoefficients_[cellEdges_[k][i]] += edgeCoeffs[i];
  }

  edgeCoefficientsUpToDate_ = true;

  return;
}
// =============================================================================
DoubleVector
StkMesh::
getEdgeCoefficientsNumerically_(
  const Teuchos::ArrayRCP<const DoubleVector> edges
  ) const
{
  int numEdges = edges.size();

  // Build an equation system for the edge coefficients alpha_k.
  // They fulfill
  //
  //    |simplex| * <u,v> = \sum_{edges e_i} alpha_i <u,e_i> <e_i,v>
  //
  // for any pair of vectors u, v in the plane of the triangle.
  //
  double vol;
  switch ( numEdges )
  {
    case 3:
      vol = getTriangleArea_( edges[0], edges[1] );
      break;
    case 6:
      // TODO Come up with a cleaner solution here.
      try
      {
        vol = getTetrahedronVolume_( edges[0], edges[1], edges[2] );
      }
      catch(...)
      {
        // If computing the volume throws an exception, then the edges
        // chosen happened to be conplanar. Changing one of those
        // fixes this.
        vol = getTetrahedronVolume_( edges[0], edges[1], edges[3] );
      }
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT_MSG( true,
                           "Can only handle triangles and tetrahedra." );
  }

  Teuchos::RCP<Teuchos::SerialSymDenseMatrix<int, double> > A =
    Teuchos::rcp( new Teuchos::SerialSymDenseMatrix<int, double>( numEdges ) );
  Teuchos::RCP<DoubleVector> rhs =
    Teuchos::rcp( new DoubleVector( numEdges ) );
  Teuchos::RCP<DoubleVector> alpha =
    Teuchos::rcp( new DoubleVector( numEdges ) );

  // Build the equation system:
  // The equation
  //
  //    |simplex| ||u||^2 = \sum_i \alpha_i <u,e_i> <e_i,u>
  //
  // has to hold for all vectors u in the plane spanned by the edges,
  // particularly by the edges themselves.
  //
  // Only fill the upper part of the Hermitian matrix.
  //
  for ( int i=0; i<numEdges; i++ )
  {
    (*rhs)(i) = vol * edges[i].dot( edges[i] );
    for ( int j=i; j<numEdges; j++ )
      (*A)(i,j) = edges[i].dot( edges[j] ) * edges[j].dot( edges[i] );
  }

  // Solve the equation system for the alpha_i.
  // The system is symmetric and, if the simplex is
  // not degenerate, positive definite.
  Teuchos::SerialSpdDenseSolver<int,double> solver;
  TEUCHOS_ASSERT_EQUALITY( 0, solver.setMatrix( A ) );
  TEUCHOS_ASSERT_EQUALITY( 0, solver.setVectors( alpha, rhs ) );
  if ( solver.shouldEquilibrate() )
  {
    TEUCHOS_ASSERT_EQUALITY( 0, solver.equilibrateRHS() );
#if TRILINOS_MAJOR_MINOR_VERSION > 100803
    TEUCHOS_ASSERT_EQUALITY( 0, solver.equilibrateMatrix() );
    TEUCHOS_ASSERT_EQUALITY( 0, solver.solve() );
    TEUCHOS_ASSERT_EQUALITY( 0, solver.unequilibrateLHS() );
#else
    // A bug in Trilinos<=10.8.3 makes unequilibrateLHS() fail.
    // Work around by doing the unequilibration manually.
    // Note that this relies on the scaling being
    // s(i) = 1 / sqrt(A(i,i)).
    Teuchos::RCP<DoubleVector> diagA =
      Teuchos::rcp( new DoubleVector( numEdges ) );
    for ( int k=0; k<numEdges; k++ )
      (*diagA)(k) = (*A)(k,k);
    TEUCHOS_ASSERT_EQUALITY( 0, solver.equilibrateMatrix() );
    TEUCHOS_ASSERT_EQUALITY( 0, solver.solve() );
    for ( int k=0; k<numEdges; k++ )
      (*alpha)(k) = (*alpha)(k) / sqrt( (*diagA)(k) );
#endif
  }
  else
  {
    TEUCHOS_ASSERT_EQUALITY( 0, solver.solve() );
  }

//     Teuchos::ArrayRCP<double> alphaArrayRCP( numEdges );
//     for ( int k=0; k<numEdges; k++ )
//         alphaArrayRCP[k] = (*alpha)(k,0);

  // TODO Don't explicitly copy alpha on return.
  return *alpha;
}
// =============================================================================
void
StkMesh::
computeControlVolumes_() const
{
  // Compute the volume of the (Voronoi) control cells for each point.
#ifndef NDEBUG
  TEUCHOS_ASSERT( !controlVolumes_.is_null() );
  TEUCHOS_ASSERT( !nodesOverlapMap_.is_null() );
#endif

  if ( !controlVolumes_->Map().SameAs( *nodesMap_ ) )
    TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->ReplaceMap( *nodesMap_ ) );
  TEUCHOS_ASSERT_EQUALITY( 0, controlVolumes_->PutScalar( 0.0 ) );

  // Create temporaries to hold the overlap values for control volumes and
  // average thickness.
  Teuchos::RCP<Epetra_Vector> cvOverlap =
    Teuchos::rcp( new Epetra_Vector( *nodesOverlapMap_ ) );

  // Determine the kind of mesh by the first cell.
  int nodesPerCell;
  if (comm_.MyPID() == 0)
  {
    std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
    nodesPerCell = cells[0]->relations(metaData_.node_rank()).size();
  }
  comm_.Broadcast(&nodesPerCell, 1, 0);

  switch (nodesPerCell)
  {
    case 3:
      this->computeControlVolumesTri_(cvOverlap);
      break;
    case 4:
      this->computeControlVolumesTet_(cvOverlap);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                                  "Illegal cell type.");
  }

  // Export control volumes to a non-overlapping map, and sum the entries.
  Epetra_Export exporter( *nodesOverlapMap_, *nodesMap_ );
  TEUCHOS_ASSERT_EQUALITY( 0,
                           controlVolumes_->Export( *cvOverlap, exporter, Add ) );

  controlVolumesUpToDate_ = true;

  return;
}
// =============================================================================
void
StkMesh::
computeControlVolumesTri_(const Teuchos::RCP<Epetra_Vector> & cvOverlap) const
{
  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  // Calculate the contributions to the finite volumes cell by cell.
  for (unsigned int k=0; k < numCells; k++)
  {
    stk::mesh::PairIterRelation localNodes =
      cells[k]->relations( metaData_.node_rank() );
    unsigned int numLocalNodes = localNodes.size();

#ifndef NDEBUG
    // Confirm that we always have the same simplices.
    TEUCHOS_ASSERT_EQUALITY(numLocalNodes, 3);
#endif

    // Fetch the nodal positions into 'localNodes'.
    //const Teuchos::ArrayRCP<const DoubleVector> localNodeCoords =
      //this->getNodeCoordinates_( localNodes );
    Teuchos::ArrayRCP<DoubleVector> localNodeCoords(numLocalNodes);
    for (unsigned int i=0; i<numLocalNodes; i++)
      localNodeCoords[i] = this->getVectorFieldNonconst(localNodes[i].entity(),
                                                        "coordinates", 3);

    // compute the circumcenter of the cell
    DoubleVector cc =
      this->computeTriangleCircumcenter_( localNodeCoords );

    // Iterate over the edges.
    // As true edge entities are not available here, loop over all pairs of local nodes.
    for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
    {
      const DoubleVector &x0 = localNodeCoords[e0];
      const int gid0 = (*localNodes[e0].entity()).identifier() - 1;
      const int lid0 = nodesOverlapMap_->LID( gid0 );
#ifndef NDEBUG
      TEUCHOS_ASSERT_INEQUALITY( lid0, >=, 0 );
#endif
      for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
      {
        const DoubleVector &x1 = localNodeCoords[e1];
        const int gid1 = (*localNodes[e1].entity()).identifier() - 1;
        const int lid1 = nodesOverlapMap_->LID( gid1 );
#ifndef NDEBUG
        TEUCHOS_ASSERT_INEQUALITY( lid1, >=, 0 );
#endif
        // Get the other node.
        Teuchos::Tuple<unsigned int,2> other = this->getOtherIndices_( e0, e1 );

        double edgeLength = this->norm2_( this->add_( 1.0, x1, -1.0, x0 ) );

        // Compute the (n-1)-dimensional covolume.
        double covolume;
        const DoubleVector &other0 = localNodeCoords[other[0]];
        covolume = this->computeCovolume2d_( cc, x0, x1, other0 );
        // The problem with counting the average thickness in 2D is the following.
        // Ideally, one would want to loop over all edges, add the midpoint value
        // of the thickness to both of the edge end points, and eventually loop over
        // all endpoints and divide by the number of edges (connections) they have
        // with neighboring nodes).
        // Unfortunately, this is impossible now b/c there's no edge generation
        // for shells in Trilinos yet (2011-04-15).
        // As a workaround, one could loop over all cells, and then all pairs of
        // nodes to retrieve the edges. In 2D, almost all of the edges would be
        // counted twice this way as they belong to two cells. This is true for
        // all but the boundary edges. Again, it is difficult (impossible?) to
        // know what the boundary edges are, and hence which values to divide by
        // 2. Dividing them all by two would result in an artificially lower
        // thickness near the boundaries. This is not what we want.

        // Compute the contributions to the finite volumes of the adjacent edges.
        double pyramidVolume = 0.5*edgeLength * covolume / 2;
        (*cvOverlap)[lid0] += pyramidVolume;
        (*cvOverlap)[lid1] += pyramidVolume;
      }
    }
  }

  return;
}
// =============================================================================
void
StkMesh::
computeControlVolumesTet_(const Teuchos::RCP<Epetra_Vector> & cvOverlap) const
{
  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numCells = cells.size();

  // Calculate the contributions to the finite volumes cell by cell.
  for (unsigned int k=0; k < numCells; k++)
  {
    stk::mesh::PairIterRelation localNodes =
      cells[k]->relations( metaData_.node_rank() );
    unsigned int numLocalNodes = localNodes.size();
#ifndef NDEBUG
    // Confirm that we always have the same simplices.
    TEUCHOS_ASSERT_EQUALITY(numLocalNodes, 4);
#endif

    // Fetch the nodal positions into 'localNodes'.
    Teuchos::ArrayRCP<DoubleVector> localNodeCoords(numLocalNodes);
    for (unsigned int i=0; i<numLocalNodes; i++)
      localNodeCoords[i] = this->getVectorFieldNonconst(localNodes[i].entity(),
                                                        "coordinates", 3);

    // compute the circumcenter of the cell
    const DoubleVector cc =
      this->computeTetrahedronCircumcenter_( localNodeCoords );

    // Iterate over the edges.
    // As true edge entities are not available here, loop over all pairs of local nodes.
    for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
    {
      const DoubleVector &x0 = localNodeCoords[e0];
      const int gid0 = (*localNodes[e0].entity()).identifier() - 1;
      const int lid0 = nodesOverlapMap_->LID( gid0 );
#ifndef NDEBUG
      TEUCHOS_ASSERT_INEQUALITY( lid0, >=, 0 );
#endif
      for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
      {
        const DoubleVector &x1 = localNodeCoords[e1];
        const int gid1 = (*localNodes[e1].entity()).identifier() - 1;
        const int lid1 = nodesOverlapMap_->LID( gid1 );
#ifndef NDEBUG
        TEUCHOS_ASSERT_INEQUALITY( lid1, >=, 0 );
#endif

        // Get the other nodes.
        Teuchos::Tuple<unsigned int,2> other = this->getOtherIndices_( e0, e1 );

        double edgeLength = this->norm2_( this->add_( 1.0, x1, -1.0, x0 ) );

        // Compute the (n-1)-dimensional covolume.
        const DoubleVector &other0 = localNodeCoords[other[0]];
        const DoubleVector &other1 = localNodeCoords[other[1]];
        double covolume = this->computeCovolume3d_( cc, x0, x1, other0, other1 );
        // Throw an exception for 3D volumes.
        // To compute the average of the thicknesses of a control volume, one has to loop
        // over all the edges and add the thickness value to both endpoints.
        // Then eventually, for each node, divide the resulting sum by the number of connections
        // (=number of faces of the finite volume).
        // However, looping over edges is not (yet) possible. Hence, we loop over all the
        // cells here. This way, the edges are counted several times, but it is difficult
        // to determine how many times exactly.
//                   TEUCHOS_TEST_FOR_EXCEPTION( true,
//                                       std::runtime_error,
//                                       "Cannot calculate the average thickness in a 3D control volume yet."
//                                     );

        // Compute the contributions to the finite volumes of the adjacent edges.
        double pyramidVolume = 0.5*edgeLength * covolume / 3;
        (*cvOverlap)[lid0] += pyramidVolume;
        (*cvOverlap)[lid1] += pyramidVolume;
      }
    }
  }
}
// =============================================================================
double
StkMesh::
computeCovolume2d_( const DoubleVector &cc,
                    const DoubleVector &x0,
                    const DoubleVector &x1,
                    const DoubleVector &other0
                    ) const
{
  // edge midpoint
  DoubleVector mp = this->add_( 0.5, x0, 0.5, x1 );

  double coedgeLength = this->norm2_( this->add_( 1.0, mp, -1.0, cc ) );

  // The only difficulty here is to determine whether the length of coedge
  // is to be taken positive or negative.
  // To this end, make sure that the order (x0, cc, mp) is of the same
  // orientation as (x0, other0, mp).
  DoubleVector cellNormal = this->cross_( this->add_( 1.0, other0, -1.0, x0 ),
                                          this->add_( 1.0, mp,     -1.0, x0 )
                                          );
  DoubleVector ccNormal = this->cross_( this->add_( 1.0, cc, -1.0, x0 ),
                                        this->add_( 1.0, mp, -1.0, x0 )
                                        );

  // copysign takes the absolute value of the first argument and the sign of the second.
  return copysign( coedgeLength, this->dot_( ccNormal, cellNormal ) );
}
// =============================================================================
double
StkMesh::
computeCovolume3d_( const DoubleVector &cc,
                    const DoubleVector &x0,
                    const DoubleVector &x1,
                    const DoubleVector &other0,
                    const DoubleVector &other1
                    ) const
{
  double covolume = 0.0;

  // edge midpoint
  DoubleVector mp = this->add_( 0.5, x0, 0.5, x1 );

  // Compute the circumcenters of the adjacent faces.
  // This could be precomputed as well.
  DoubleVector ccFace0 = this->computeTriangleCircumcenter_( x0, x1, other0 );
  DoubleVector ccFace1 = this->computeTriangleCircumcenter_( x0, x1, other1 );

  // Compute the area of the quadrilateral.
  // There are some really tricky degenerate cases here, i.e., combinations
  // of when ccFace{0,1}, cc, sit outside of the tetrahedron.

  // Use the triangle (MP, localNodes[other[0]], localNodes[other[1]] ) (in this order)
  // to gauge the orientation of the two triangles that compose the quadrilateral.
  DoubleVector gauge = this->cross_( this->add_( 1.0, other0, -1.0, mp ),
                                     this->add_( 1.0, other1, -1.0, mp )
                                     );

  // Add the area of the first triangle (MP,ccFace0,cc).
  // This makes use of the right angles.
  double triangleHeight0 = this->norm2_( this->add_( 1.0, mp, -1.0, ccFace0 ) );
  double triangleArea0 = 0.5
    * triangleHeight0
    * this->norm2_( this->add_( 1.0, ccFace0, -1.0, cc ) );

  // Check if the orientation of the triangle (MP,ccFace0,cc) coincides with
  // the orientation of the gauge triangle. If yes, add the area, subtract otherwise.
  DoubleVector triangleNormal0 =
    this->cross_(this->add_(1.0, ccFace0, -1.0, mp),
                 this->add_(1.0, cc,      -1.0, mp));

  // copysign takes the absolute value of the first argument and the sign of the second.
  covolume += copysign( triangleArea0, this->dot_( triangleNormal0, gauge ) );

  // Add the area of the second triangle (MP,cc,ccFace1).
  // This makes use of the right angles.
  double triangleHeight1 = this->norm2_( this->add_( 1.0, mp, -1.0, ccFace1 ) );
  double triangleArea1 = 0.5
    * triangleHeight1
    * this->norm2_( this->add_( 1.0, ccFace1, -1.0, cc ) );

  // Check if the orientation of the triangle (MP,cc,ccFace1) coincides with
  // the orientation of the gauge triangle. If yes, add the area, subtract otherwise.
  DoubleVector triangleNormal1 =
    this->cross_( this->add_(1.0, cc,      -1.0, mp),
                  this->add_(1.0, ccFace1, -1.0, mp));

  // copysign takes the absolute value of the first argument and the sign of the second.
  covolume += copysign( triangleArea1, this->dot_( triangleNormal1, gauge ) );

  return covolume;
}
// =============================================================================
Teuchos::Tuple<unsigned int,2>
StkMesh::
getOtherIndices_( unsigned int e0, unsigned int e1 ) const
{
  // Get the two indices in [0,1,2,3] which are not e0, e1.
  unsigned int count = 0;
  Teuchos::Tuple<unsigned int,2> otherInd;
  for ( unsigned int k=0; k<4; k++ )
  {
    if ( k!=e0 && k!=e1 )
      otherInd[count++] = k;
#ifndef NDEBUG
    TEUCHOS_ASSERT_INEQUALITY( count, <=, 2 );
#endif
  }
  return otherInd;
}
// =============================================================================
double
StkMesh::
getTriangleArea_( const DoubleVector &edge0,
                  const DoubleVector &edge1
                  ) const
{
  return 0.5 * this->norm2_( this->cross_( edge0, edge1 ) );
}
// =============================================================================
double
StkMesh::
getTetrahedronVolume_( const DoubleVector &edge0,
                       const DoubleVector &edge1,
                       const DoubleVector &edge2
                       ) const
{
  // Make sure the edges are not conplanar.
  double alpha = edge0.dot( this->cross_( edge1, edge2 ) );
  TEUCHOS_TEST_FOR_EXCEPT_MSG( fabs( alpha ) / this->norm2_( edge0 )
                       / this->norm2_( edge1 )
                       / this->norm2_( edge2 )
                       < 1.0e-5,
                       "The following edges seem to be conplanar:"
                       << "\n(0) " << edge0
                       << "\n(1) " << edge1
                       << "\n(2) " << edge2 );
  double vol = fabs( alpha ) / 6.0;
  return vol;
}
// =============================================================================
DoubleVector
StkMesh::
computeTriangleCircumcenter_(
  const Teuchos::ArrayRCP<const DoubleVector> &nodes ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY( nodes.size(), 3 );
#endif
  return this->computeTriangleCircumcenter_( nodes[0], nodes[1], nodes[2] );
}
// =============================================================================
DoubleVector
StkMesh::
computeTriangleCircumcenter_( const DoubleVector &node0,
                              const DoubleVector &node1,
                              const DoubleVector &node2
                              ) const
{
  DoubleVector a(node0);
  for (int k=0; k<3; k++)
    a[k] -= node1[k];

  DoubleVector b(node1);
  for (int k=0; k<3; k++)
    b[k] -= node2[k];

  DoubleVector c(node2);
  for (int k=0; k<3; k++)
    c[k] -= node0[k];

  const double omega = 2.0 * this->norm2squared_(this->cross_(a, b));

  // don't divide by 0
  TEUCHOS_TEST_FOR_EXCEPT_MSG( fabs( omega ) < 1.0e-10,
                       "It seems that the nodes \n\n"
                       << "   " << node0 << "\n"
                       << "   " << node1 << "\n"
                       << "   " << node2 << "\n"
                       << "\ndo not form a proper triangle. Abort."
                       << std::endl
                       );

  const double alpha = - this->dot_(b, b) * this->dot_(a, c) / omega;
  const double beta  = - this->dot_(c, c) * this->dot_(b, a) / omega;
  const double gamma = - this->dot_(a, a) * this->dot_(c, b) / omega;

  DoubleVector cc(3);
  for (int k=0; k<3; k++)
    cc[k] = alpha * node0[k]
          + beta  * node1[k]
          + gamma * node2[k];

  return cc;
}
// =============================================================================
DoubleVector
StkMesh::
computeTetrahedronCircumcenter_(
  const Teuchos::ArrayRCP<const DoubleVector> &nodes ) const
{
  // http://www.cgafaq.info/wiki/Tetrahedron_Circumsphere
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY( nodes.size(), 4 );
#endif

  // Compute with respect to the first point.
  Teuchos::Array<DoubleVector> relNodes( 3 );
  for ( int k=0; k<3; k++ )
    relNodes[k] = this->add_( 1.0, nodes[k+1], -1.0, nodes[0] );


  double omega = 2.0 *
    this->dot_( relNodes[0], this->cross_( relNodes[1], relNodes[2] ) );

  // don't divide by 0
  TEUCHOS_TEST_FOR_EXCEPT_MSG( fabs( omega ) < 1.0e-10,
                       "It seems that the nodes \n\n"
                       << "   " << nodes[0] << "\n"
                       << "   " << nodes[1] << "\n"
                       << "   " << nodes[2] << "\n"
                       << "   " << nodes[3] << "\n"
                       << "\ndo not form a proper tetrahedron. Abort."
                       << std::endl
                       );
  double alpha = this->norm2squared_( relNodes[0] ) / omega;
  double beta  = this->norm2squared_( relNodes[1] ) / omega;
  double gamma = this->norm2squared_( relNodes[2] ) / omega;

  DoubleVector cc;
  cc = this->add_( alpha, this->cross_( relNodes[1], relNodes[2] ),
                   beta,  this->cross_( relNodes[2], relNodes[0] ) );
  cc = this->add_( 1.0,   cc,
                   gamma, this->cross_( relNodes[0], relNodes[1] ) );

  cc = this->add_( 1.0, cc,
                   1.0, nodes[0] );

  return cc;
}
// =============================================================================
DoubleVector
StkMesh::
add_( double alpha, const DoubleVector &x,
      double beta,  const DoubleVector &y
      ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY( x.length(), 3 );
  TEUCHOS_ASSERT_EQUALITY( y.length(), 3 );
#endif
  DoubleVector z( 3 );
  for ( int k=0; k<z.length(); k++ )
    z[k] = alpha*x[k] + beta*y[k];

  return z;
}
// =============================================================================
double
StkMesh::
dot_( const DoubleVector &v,
      const DoubleVector &w
      ) const
{
  double sum = 0.0;
  for ( int k=0; k<v.length(); k++ )
    sum += v[k] * w[k];
  return sum;
}
// =============================================================================
DoubleVector
StkMesh::
cross_( const DoubleVector &v,
        const DoubleVector &w
        ) const
{
  DoubleVector z( 3 );

  z[0] = v[1]*w[2] - v[2]*w[1];
  z[1] = v[2]*w[0] - v[0]*w[2];
  z[2] = v[0]*w[1] - v[1]*w[0];

  return z;
}
// =============================================================================
double
StkMesh::
norm2_( const DoubleVector &x
        ) const
{
  return sqrt( this->dot_( x,x ) );
}
// =============================================================================
double
StkMesh::
norm2squared_( const DoubleVector &x
               ) const
{
  return this->dot_( x, x );
}
// =============================================================================
void
StkMesh::
createEdges_()
{
  if (edgeNodes_.size() != 0 && !cellEdges_.is_null())
    return;

  std::vector<stk::mesh::Entity*> cells = this->getOwnedCells();
  unsigned int numLocalCells = cells.size();
  // Local cell ID -> Local edge IDs.
  cellEdges_ = Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> >(numLocalCells);

  // This std::map keeps track of how nodes and edges are connected.
  // If  nodeEdges((3,4)) == 17  is true, then the nodes (3,4) are
  // connected  by edge 17.
  // Unfortunately, Teuchos::Tuples can't be compared with '<'. Provide a
  // function pointer that implements lexicographic comparison.
  // See http://www.cplusplus.com/reference/stl/map/map/.
  std::map<Teuchos::Tuple<stk::mesh::Entity*,2>,int,TupleComp> nodesEdge;

  // Loop over all owned cells.
  unsigned int edgeLID = 0;
  for ( unsigned int cellLID=0; cellLID<numLocalCells; cellLID++ )
  {
    // Loop over all pairs of local nodes.
    stk::mesh::PairIterRelation nodesIterator =
      cells[cellLID]->relations( metaData_.node_rank() );
    unsigned int numLocalNodes = nodesIterator.size();
    unsigned int numLocalEdges = numLocalNodes*(numLocalNodes-1) / 2;

    cellEdges_[cellLID] = Teuchos::ArrayRCP<int>(numLocalEdges);

    // Gather the node entities.
    Teuchos::ArrayRCP<stk::mesh::Entity*> nodes(numLocalNodes);
    for ( unsigned int k=0; k<numLocalNodes; k++ )
      nodes[k] = nodesIterator[k].entity();

    // Sort nodes by their global identifier. This is necessary
    // to make sure that the tuples formed below are always sorted
    // such they are unique keys (and {3,7}, {7,3} are recognized
    // as the same edge).
    EntityComp ec;
    std::sort(nodes.begin(), nodes.end(), ec);

    // In a simplex, the edges are exactly the connection between each pair
    // of nodes. Hence, loop over pairs of nodes.
    unsigned int edgeIndex = 0;
    Teuchos::Tuple<stk::mesh::Entity*,2> edgeNodes;
    for ( unsigned int e0=0; e0<numLocalNodes; e0++ )
    {
      edgeNodes[0] = nodes[e0];
      for ( unsigned int e1=e0+1; e1<numLocalNodes; e1++ )
      {
        edgeNodes[1] = nodes[e1];
        // As nodes are sorted and by their identifiers, edgeNodes are sorted
        // too. This is necessary as otherwise the edge {3,7} could not be
        // identified as {7,3}.

        // Check if edgeNodes is in the map.
        std::map<Teuchos::Tuple<stk::mesh::Entity*,2>,int,TupleComp>::iterator it =
            nodesEdge.find(edgeNodes);
        if ( it != nodesEdge.end() )
        {
          // Edge is already accounted for.
          cellEdges_[cellLID][edgeIndex] = it->second;
        }
        else // Edge not found -- insert it.
        {
          nodesEdge[edgeNodes] = edgeLID; // for householding in this method
          edgeNodes_.append( edgeNodes ); // for looping over edges
          cellEdges_[cellLID][edgeIndex] = edgeLID; // for this->computeEdgeCoefficients_
          edgeLID++;
        }
        edgeIndex++;
      }
    }
  }

  return;
}
// =============================================================================
} // namespace Nosh
