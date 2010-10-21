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

#include "Ginla_EpetraFVM_State.h"

#include "VIO_EpetraMesh_Writer.h"

#include <LOCA_Parameter_Vector.H>
#include <Epetra_Comm.h>

// =============================================================================
Ginla::EpetraFVM::State::
State( const Epetra_Vector                             & psi,
       const Teuchos::RCP<const VIO::EpetraMesh::Mesh> & mesh
     ):
       psi_( psi ),
       mesh_( mesh )
{
    TEUCHOS_ASSERT_EQUALITY( psi.GlobalLength(),
                             2 * mesh->getNumNodes() );
}
// =============================================================================
Ginla::EpetraFVM::State::
State( const Teuchos::RCP<const Epetra_Map>            & map,
       const Teuchos::RCP<const VIO::EpetraMesh::Mesh> & mesh
     ):
       psi_( Epetra_Vector( *map, true ) ),
       mesh_( mesh )
{
}
// =============================================================================
Ginla::EpetraFVM::State::
State( const Teuchos::RCP<const Epetra_Comm>           & comm,
       const Teuchos::RCP<const VIO::EpetraMesh::Mesh> & mesh
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
const Teuchos::RCP<const VIO::EpetraMesh::Mesh>
Ginla::EpetraFVM::State::
getMesh () const
{
    return mesh_;
}
// =============================================================================
void
Ginla::EpetraFVM::State::
save( const std::string            & fileName,
      const Teuchos::ParameterList & p
    ) const
{
    TEUCHOS_ASSERT( !mesh_.is_null() );

    Teuchos::RCP<VIO::EpetraMesh::Writer> writer =
        Teuchos::rcp( new VIO::EpetraMesh::Writer( fileName ) );

    writer->setMesh( *mesh_ );
    writer->setValues( psi_ );
    writer->addParameterList( p );

    writer->write();

    return;
}
// =============================================================================
void
Ginla::EpetraFVM::State::
save( const std::string & fileName
    ) const
{
    Teuchos::ParameterList empty;
    this->save( fileName, empty );
}
// =============================================================================
double
Ginla::EpetraFVM::State::
freeEnergy () const
{
    double myGlobalEnergy[1];

    int numMyPoints = mesh_->getNodesMap()->NumMyPoints();
    TEUCHOS_ASSERT_EQUALITY( 2*numMyPoints, psi_.MyLength() );

    const Epetra_Vector & controlVolumes =  *(mesh_->getControlVolumes());

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

    int numMyPoints = mesh_->getNodesMap()->NumMyPoints();
    TEUCHOS_ASSERT_EQUALITY( 2*numMyPoints, psi_.MyLength() );

    const Epetra_Vector & controlVolumes = *(mesh_->getControlVolumes());
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
// ============================================================================
