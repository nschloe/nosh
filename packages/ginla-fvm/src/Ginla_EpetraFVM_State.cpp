/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Schl"omer

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
#include "Ginla_EpetraFVM_State.h"

#include "Ginla_EpetraFVM_StkMesh.h"
#include "Ginla_EpetraFVM_StkMeshWriter.h"

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <LOCA_Parameter_Vector.H>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
// =============================================================================
Ginla::EpetraFVM::State::
State( const Epetra_Vector                                 & psi,
       const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh
     ):
       psi_( psi ),
       mesh_( mesh )
{
    TEUCHOS_ASSERT_EQUALITY( psi.GlobalLength(),
                             2 * mesh->getNumNodes() );

    return;
}
// =============================================================================
Ginla::EpetraFVM::State::
State( const Teuchos::RCP<const Epetra_Map>                & map,
       const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh
     ):
       psi_( Epetra_Vector( *map, true ) ),
       mesh_( mesh )
{
}
// =============================================================================
Ginla::EpetraFVM::State::
State( const Teuchos::RCP<const Epetra_Comm>               & comm,
       const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh
     ):
       psi_( Epetra_Vector( Epetra_Map( 2*mesh->getNumNodes(), 0, *comm ) ) ),
       mesh_( mesh )
{
}
// =============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::EpetraFVM::State::
getPsi () const
{
    return Teuchos::rcpFromRef( psi_ );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::EpetraFVM::State::
getPsiNonConst ()
{
return Teuchos::rcpFromRef( psi_ );
}
// =============================================================================
const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh>
Ginla::EpetraFVM::State::
getMesh () const
{
    return mesh_;
}
// =============================================================================
void
Ginla::EpetraFVM::State::
save( const std::string            & fileBaseName,
      const int                      index,
      const Teuchos::ParameterList & p
    ) const
{

    // Merge the state into the mesh.
//     mesh_->getBulkData()->modification_begin();
    this->mergePsi_( mesh_, psi_ );
//     mesh_->getBulkData()->modification_end();


    // Write it out to the file that's been specified previously.
    int out_step = stk::io::util::process_output_request( *mesh_->getMeshData(),
                                                          *mesh_->getBulkData(),
                                                          index
                                                        );

    std::cout << "Ginla::EpetraFVM::StkMeshWriter::write:\n"
              << "\twriting time " << index << "\n"
              << "\tindex " << out_step << "\n"
              << "\tto file " //<< fileName.str()
              << std::endl;

//     Ginla::EpetraFVM::StkMeshWrite( fileBaseName, index,  psi_, mesh_, p );

    return;
}
// =============================================================================
void
Ginla::EpetraFVM::State::
save( const std::string            & fileBaseName,
      const int                      index
    ) const
{
    Teuchos::ParameterList empty;
    this->save( fileBaseName, index, empty );
}
// =============================================================================
double
Ginla::EpetraFVM::State::
freeEnergy () const
{
    double myGlobalEnergy[1];

    const Epetra_Vector & controlVolumes =  *(mesh_->getControlVolumes());

    int numMyPoints = controlVolumes.Map().NumMyPoints();
    TEUCHOS_ASSERT_EQUALITY( 2*numMyPoints, psi_.MyLength() );

    for ( int k=0; k<numMyPoints; k++ )
        myGlobalEnergy[0] -= controlVolumes[k] * pow( psi_[2*k]*psi_[2*k] + psi_[2*k+1]*psi_[2*k+1], 2 );

    // Sum over all processors.
    const Epetra_Comm & comm = psi_.Comm();
    double globalEnergy[1];
    TEUCHOS_ASSERT_EQUALITY( 0, comm.SumAll( myGlobalEnergy, globalEnergy, 1 ) );

    globalEnergy[0] /= mesh_->getDomainArea();

    return globalEnergy[0];
}
// =============================================================================
double
Ginla::EpetraFVM::State::
innerProduct( const Ginla::EpetraFVM::State & state ) const
{
    double res[1];

    const Epetra_Vector & controlVolumes = *(mesh_->getControlVolumes());

    int numMyPoints = controlVolumes.Map().NumMyPoints();
    TEUCHOS_ASSERT_EQUALITY( 2*numMyPoints, psi_.MyLength() );

    const Epetra_Vector & psi2 =  *(state.getPsi());

    for ( int k=0; k<numMyPoints; k++ )
        res[0] += controlVolumes[k] * ( psi_[2*k]*psi2[2*k] + psi_[2*k+1]*psi2[2*k+1]  );

    // Sum over all processors.
    const Epetra_Comm & comm = psi_.Comm();
    double globalRes[1];
    TEUCHOS_ASSERT_EQUALITY( 0, comm.SumAll( res, globalRes, 1 ) );

    // normalize and return
    return globalRes[0] / mesh_->getDomainArea();
}
//// ============================================================================
double
Ginla::EpetraFVM::State::
normalizedScaledL2Norm () const
{
    return sqrt ( this->innerProduct( *this ) );
}
// =============================================================================
void
Ginla::EpetraFVM::State::
update( const double                    alpha,
        const Ginla::EpetraFVM::State & b,
        const double                    beta
      )
{
  psi_.Update( alpha, *(b.getPsi()), beta );
  return;
}
// =============================================================================
// Count the number of vortices by the total phase change along the boundary
// of the domain.
// TODO Make this work in multicore environments.
// Idea: Cauchy's integral formula: just caluculate the integral.
//       Numerically difficult when too close to origin (rather: points
//       closest to and furthest from origin too far apart).
int
Ginla::EpetraFVM::State::
getVorticity () const
{
//   TEST_FOR_EXCEPTION( true,
//                       std::logic_error,
//                       "Method \"getVorticity()\" not yet implemented."
//                     );

  return 0;
}
// =============================================================================
void
Ginla::EpetraFVM::State::
mergePsi_( const Teuchos::RCP<const Ginla::EpetraFVM::StkMesh> & mesh,
           const Epetra_Vector                                 & psi
         ) const
{
    // Get owned nodes.
    const std::vector<stk::mesh::Entity*> & ownedNodes =  mesh->getOwnedNodes();

    VectorFieldType * psi_field = mesh->getMetaData()->get_field<VectorFieldType>( "psi" );
    TEUCHOS_ASSERT( psi_field != 0 );

    // Merge psi into the mesh.
    for (int k=0; k < ownedNodes.size(); k++)
    {
        // Extract real and imaginary part.
        double* localPsi = stk::mesh::field_data( *psi_field, *ownedNodes[k] );
        localPsi[0] = psi[2*k];
        localPsi[1] = psi[2*k+1];
    }

    // This communication updates the field values on un-owned nodes
    // it is correct because the zeroSolutionField above zeros them all
    // and the getSolutionField only sets the owned nodes.
    stk::mesh::parallel_reduce( *mesh->getBulkData() , stk::mesh::sum(*psi_field));

    return;
}
// ============================================================================
