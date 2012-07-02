// @HEADER
//
//    Container class for quantum states.
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
#include "Cuantico_State.hpp"

#include "Cuantico_StkMesh.hpp"

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>

#include <LOCA_Parameter_Vector.H>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

#ifdef CUANTICO_TEUCHOS_TIME_MONITOR
  #include <Teuchos_TimeMonitor.hpp>
#endif

namespace Cuantico {
// =============================================================================
State::
State( const Epetra_Vector &psi,
       const Teuchos::RCP<const Cuantico::StkMesh> &mesh
       ) :
#ifdef CUANTICO_TEUCHOS_TIME_MONITOR
  saveTime_( Teuchos::TimeMonitor::getNewTimer( "Cuantico: State::save" ) ),
#endif
  psi_( psi ),
  mesh_( mesh )
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( psi.GlobalLength(),
                           2 * mesh->getNumNodes() );
#endif

  return;
}
// =============================================================================
State::
State( const Teuchos::RCP<const Epetra_Map> &map,
       const Teuchos::RCP<const Cuantico::StkMesh> &mesh
       ) :
  psi_( Epetra_Vector( *map, true ) ),
  mesh_( mesh )
{
}
// =============================================================================
State::
State( const Teuchos::RCP<const Epetra_Comm> &comm,
       const Teuchos::RCP<const Cuantico::StkMesh> &mesh
       ) :
  psi_( Epetra_Vector( Epetra_Map( 2*mesh->getNumNodes(), 0, *comm ) ) ),
  mesh_( mesh )
{
}
// =============================================================================
Teuchos::RCP<const Epetra_Vector>
State::
getPsi() const
{
  return Teuchos::rcpFromRef( psi_ );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
State::
getPsiNonConst()
{
  return Teuchos::rcpFromRef( psi_ );
}
// =============================================================================
const Teuchos::RCP<const Cuantico::StkMesh>
State::
getMesh() const
{
  return mesh_;
}
// =============================================================================
void
State::
save( const int index
      ) const
{
#ifdef CUANTICO_TEUCHOS_TIME_MONITOR
  // timer for this routine
  Teuchos::TimeMonitor tm( *saveTime_ );
#endif

  // Merge the state into the mesh.
//     mesh_->getBulkData()->modification_begin();
  this->mergePsi_(mesh_, psi_);
//     mesh_->getBulkData()->modification_end();

  // Write it out to the file that's been specified in mesh_.
  double time = index;
  int out_step = stk::io::process_output_request( *mesh_->getMeshData(),
                                                  *mesh_->getBulkData(),
                                                  time
                                                  );
  return;
}
// =============================================================================
double
State::
freeEnergy() const
{
  double myEnergy = 0.0;

  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());

  int numMyPoints = controlVolumes.Map().NumMyPoints();
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( 2*numMyPoints, psi_.MyLength() );
#endif

  double alpha;
  for ( int k=0; k<numMyPoints; k++ )
  {
    alpha = psi_[2*k]*psi_[2*k] + psi_[2*k+1]*psi_[2*k+1];
    myEnergy -= controlVolumes[k] * alpha * alpha;
  }

  // Sum over all processors.
  const Epetra_Comm &comm = psi_.Comm();
  double globalEnergy;
  TEUCHOS_ASSERT_EQUALITY( 0, comm.SumAll( &myEnergy, &globalEnergy, 1 ) );

  // normalize and return
  return globalEnergy / mesh_->getDomainVolume();
}
// =============================================================================
double
State::
innerProduct( const Cuantico::State &state ) const
{
  double res = 0.0;

  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());

  int numMyPoints = controlVolumes.Map().NumMyPoints();
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( 2*numMyPoints, psi_.MyLength() );
#endif

  const Epetra_Vector &psi2 = *(state.getPsi());

  for ( int k=0; k<numMyPoints; k++ )
    res += controlVolumes[k] * (psi_[2*k]*psi2[2*k] + psi_[2*k+1]*psi2[2*k+1]);

  // Sum over all processors.
  const Epetra_Comm &comm = psi_.Comm();
  double globalRes;
  TEUCHOS_ASSERT_EQUALITY( 0, comm.SumAll( &res, &globalRes, 1 ) );

  // normalize and return
  return globalRes / mesh_->getDomainVolume();
}
// ============================================================================
double
State::
normalizedScaledL2Norm() const
{
  return sqrt( this->innerProduct( *this ) );
}
// =============================================================================
void
State::
mergePsi_( const Teuchos::RCP<const Cuantico::StkMesh> &mesh,
           const Epetra_Vector &psi
           ) const
{
  VectorFieldType * psir_field =
    mesh->getMetaData()->get_field<VectorFieldType>("psi_R");
#ifdef _DEBUG_
  TEUCHOS_ASSERT( psir_field != NULL );
#endif

  VectorFieldType * psii_field =
    mesh->getMetaData()->get_field<VectorFieldType>("psi_Z");
#ifdef _DEBUG_
  TEUCHOS_ASSERT( psii_field != NULL );
#endif

  // Zero out all nodal values.
  const std::vector<stk::mesh::Entity*> &overlapNodes = mesh->getOverlapNodes();
  for (unsigned int k=0; k < overlapNodes.size(); k++)
  {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data( *psir_field, *overlapNodes[k] );
    localPsiR[0] = 0.0;
    double* localPsiI = stk::mesh::field_data( *psii_field, *overlapNodes[k] );
    localPsiI[0] = 0.0;
  }

  // Set owned nodes.
  const std::vector<stk::mesh::Entity*> &ownedNodes = mesh->getOwnedNodes();
  for (unsigned int k=0; k < ownedNodes.size(); k++)
  {
    // Extract real and imaginary part.
    double* localPsiR = stk::mesh::field_data( *psir_field, *ownedNodes[k] );
    localPsiR[0] = psi[2*k];
    double* localPsiI = stk::mesh::field_data( *psii_field, *ownedNodes[k] );
    localPsiI[0] = psi[2*k+1];
  }

  // This communication updates the field values on un-owned nodes
  // it is correct because the zeroSolutionField above zeros them all
  // and the getSolutionField only sets the owned nodes.
  stk::mesh::parallel_reduce(*mesh->getBulkData(),
                             stk::mesh::sum(*psir_field));
  stk::mesh::parallel_reduce(*mesh->getBulkData(),
                             stk::mesh::sum(*psii_field));

  return;
}
// ============================================================================
} // namespace Cuantico
