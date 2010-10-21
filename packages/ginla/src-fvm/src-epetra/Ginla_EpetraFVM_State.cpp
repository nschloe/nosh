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
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Not yet implemented."
                      );

//   double myGlobalEnergy = 0.0;
//   Teuchos::ArrayRCP<const double> cvView =
//       mesh_->getControlVolumes()->get1dView();
//   Teuchos::ArrayRCP<const double_complex> psiView = psi_.get1dView();
//   for ( int k=0; k<mesh_->getControlVolumes()-> psiView.size(); k++ )
//   {
//       // get the local index for the control volumes
//       ORD globalIndex = psi_.getMap()->getGlobalElement( k );
//       TEUCHOS_ASSERT( globalIndex != Teuchos::OrdinalTraits<ORD>::invalid() );
//       ORD k2 = mesh_->getControlVolumes()->getMap()->getLocalElement( globalIndex );
//       TEUCHOS_ASSERT( k2 != Teuchos::OrdinalTraits<ORD>::invalid() );
//
//       myGlobalEnergy -= cvView[k2] * pow( norm( psiView[k] ), 2 );
//   }
//
//   // sum over all processes
//   int count = 1; // send *one* integer
//   Teuchos::Array<double> sendBuff ( count );
//   sendBuff[0] = myGlobalEnergy;
//   Teuchos::Array<double> recvBuff ( count );
//   Teuchos::reduceAll ( * ( psi_.Map()->getComm() ),
//                        Teuchos::REDUCE_SUM,
//                        count,
//                        sendBuff.getRawPtr(),
//                        recvBuff.getRawPtr()
//                      );
//
//   recvBuff[0] /= mesh_->getDomainArea();
//
//    return recvBuff[0];

    return 0.0;
}
// =============================================================================
double_complex
Ginla::EpetraFVM::State::
innerProduct( const Ginla::EpetraFVM::State & state ) const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Not yet implemented."
                      );

//    Teuchos::ArrayRCP<const double_complex> psiView = psi_.get1dView();
//    Teuchos::ArrayRCP<const double_complex> phiView = state.getPsi()->get1dView();
//
//    // make sure the maps are the same
//    TEUCHOS_ASSERT( state.getPsi()->Map()->isSameAs( *psi_.getMap() ) );
//
//    // TODO replace by Tpetra::weighted inner product when available
//    double_complex localSum = double_complex( 0.0, 0.0 );
//    Teuchos::ArrayRCP<const double> controlVolumes =
//        mesh_->getControlVolumes()->get1dView();
//    for ( int k=0; k<psiView.size(); k++ )
//    {
//        int globalIndex = psi_.Map()->getGlobalElement( k );
//        // translate into controlVolumes local index
//        int k2 = mesh_->getControlVolumes()->getMap()->getLocalElement( globalIndex );
//        localSum += controlVolumes[k2] * conj(psiView[k]) * phiView[k];
//    }
//
//    // reduce and scatter such that energy is available on
//    // all cores
//    int count = 1; // send *one* integer
//
//    Teuchos::Array<double_complex> sendBuff ( count );
//    sendBuff[0] = localSum;
//
//    Teuchos::Array<double_complex> recvBuff ( count );
//
//    Teuchos::reduceAll ( * ( psi_.getMap()->getComm() ),
//                         Teuchos::REDUCE_SUM,
//                         count,
//                         sendBuff.getRawPtr(),
//                         recvBuff.getRawPtr()
//                       );
//
//    return recvBuff[0];
    return 0.0;
}
//// ============================================================================
double
Ginla::EpetraFVM::State::
normalizedScaledL2Norm () const
{
    // TODO replace by Tpetra::weighted norm as soon as available

    // imaginary part of alpha should be 0
    double_complex alpha = this->innerProduct( *this );

    // make sure that we actually got a norm here
    TEUCHOS_ASSERT_INEQUALITY( alpha.imag(), <, 1.0e-10 );

    // normalize
    double domainArea = mesh_->getDomainArea();
    double l2norm = sqrt ( alpha.real() ) / domainArea;

    return l2norm;
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
