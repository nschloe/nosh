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

#include "Ginla_FVM_State.h"

#include "VIO_Mesh_Writer.h"

#include <LOCA_Parameter_Vector.H>

// =============================================================================
Ginla::FVM::State::
State( const Teuchos::RCP<ComplexMultiVector>    & psi,
       const Teuchos::RCP<const VIO::Mesh::Mesh> & mesh ):
       psi_( *psi ),
       chi_( 0.0 ),
       mesh_( mesh )
{
    TEUCHOS_ASSERT_EQUALITY( psi->getGlobalLength(),
                             mesh->getNumNodes() );
}
// =============================================================================
Ginla::FVM::State::
State( const Teuchos::RCP<ComplexVector>         & psi,
       const Teuchos::RCP<const VIO::Mesh::Mesh> & mesh ):
       psi_( *psi ),
       chi_( 0.0 ),
       mesh_( mesh )
{
    TEUCHOS_ASSERT_EQUALITY( psi->getGlobalLength(),
                             mesh->getNumNodes() );
}
// =============================================================================
Ginla::FVM::State::
State( const Teuchos::RCP<const ComplexMap>      & map,
       const Teuchos::RCP<const VIO::Mesh::Mesh> & mesh ):
       psi_( ComplexMultiVector( map, 1, true ) ),
       chi_( 0.0 ),
       mesh_( mesh )
{
}
// =============================================================================
Ginla::FVM::State::
State( const Teuchos::RCP<const Teuchos::Comm<int> >  & comm,
       const Teuchos::RCP<const VIO::Mesh::Mesh>      & mesh
     ):
       psi_( ComplexMultiVector( Teuchos::rcp( new ComplexMap( mesh->getNumNodes(), 0, comm ) ),
                                 1,
                                 true ) ),
       chi_( 0.0 ),
       mesh_( mesh )
{  
}
// =============================================================================
Teuchos::RCP<const ComplexVector>
Ginla::FVM::State::
getPsi () const
{
    return psi_.getVector(0); 
}
// =============================================================================
Teuchos::RCP<ComplexVector>
Ginla::FVM::State::
getPsiNonConst ()
{
    return psi_.getVectorNonConst(0); 
}
// =============================================================================
const Teuchos::RCP<const VIO::Mesh::Mesh>
Ginla::FVM::State::
getMesh () const
{
    return mesh_; 
}
// =============================================================================
double
Ginla::FVM::State::
getChi () const
{
  return chi_;
}
// =============================================================================
void
Ginla::FVM::State::
save( const std::string            & fileName,
      const Teuchos::ParameterList & p
    ) const
{
    TEUCHOS_ASSERT( !mesh_.is_null() );
  
    Teuchos::RCP<VIO::Mesh::Writer> writer =
        Teuchos::rcp( new VIO::Mesh::Writer( fileName ) );
        
    writer->setMesh( *mesh_ );
    writer->setValues( psi_ );
    // TODO set Teuchos::ParameterList
    
    writer->write();
    
    return;
}
// =============================================================================
void
Ginla::FVM::State::
save( const std::string & fileName
    ) const
{
    Teuchos::ParameterList empty;
    this->save( fileName, empty );
}
// =============================================================================
double
Ginla::FVM::State::
freeEnergy () const
{
   double myGlobalEnergy = 0.0;
   Teuchos::ArrayRCP<const double> cvView =
       mesh_->getControlVolumes()->get1dView();
   Teuchos::ArrayRCP<const double_complex> psiView = psi_.get1dView();
   for ( int k=0; k<psiView.size(); k++ )
   {
       // get the local index for the control volumes
       ORD globalIndex = psi_.getMap()->getGlobalElement( k );
       TEUCHOS_ASSERT( globalIndex != Teuchos::OrdinalTraits<ORD>::invalid() );
       ORD k2 = mesh_->getControlVolumes()->getMap()->getLocalElement( globalIndex );
       TEUCHOS_ASSERT( k2 != Teuchos::OrdinalTraits<ORD>::invalid() );
      
       myGlobalEnergy -= cvView[k2] * pow( norm( psiView[k] ), 2 );
   }

   // sum over all processes
   int count = 1; // send *one* integer
   Teuchos::Array<double> sendBuff ( count );
   sendBuff[0] = myGlobalEnergy;
   Teuchos::Array<double> recvBuff ( count );
   Teuchos::reduceAll ( * ( psi_.getMap()->getComm() ),
                        Teuchos::REDUCE_SUM,
                        count,
                        sendBuff.getRawPtr(),
                        recvBuff.getRawPtr()
                      );
   
   recvBuff[0] /= mesh_->getDomainArea();
   

//    // FEM calculations
//    double energy = 0.0;
//    Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> > & elems = mesh_->getElems();
//    for ( int k=0; k<elems.size(); k++ )
//    {
//        Teuchos::ArrayRCP<int> & elem = elems[k];
//        
//        TEST_FOR_EXCEPTION( elem.size() != 3,
//                            std::runtime_error,
//                            "At the moment, the free energy of a state "
//                            << "can only be computed for triangular elements."
//                          );
//         // Map the linear function to the reference element, and evaluate the
//         // integral there. These calculations can be exectuted analytically,
//         // and the result is
//         //    \int_{element} psi     = |detJ| * 1/6 * ( y0 + y1 + y2 )
//         //    \int_{element} |psi|^4 = ?
// 
//    }

    return recvBuff[0];
}
// =============================================================================
double_complex
Ginla::FVM::State::
innerProduct( const Ginla::FVM::State & state ) const
{ 
    Teuchos::ArrayRCP<const double_complex> psiView = psi_.get1dView();
    Teuchos::ArrayRCP<const double_complex> phiView = state.getPsi()->get1dView();
    
    // make sure the maps are the same
    TEUCHOS_ASSERT( state.getPsi()->getMap()->isSameAs( *psi_.getMap() ) );
  
    // TODO replace by Tpetra::weighted inner product when available
    double_complex localSum = double_complex( 0.0, 0.0 );
    Teuchos::ArrayRCP<const double> controlVolumes =
        mesh_->getControlVolumes()->get1dView();
    for ( int k=0; k<psiView.size(); k++ )
    {
        int globalIndex = psi_.getMap()->getGlobalElement( k );
        // translate into controlVolumes local index
        int k2 = mesh_->getControlVolumes()->getMap()->getLocalElement( globalIndex );
        localSum += controlVolumes[k2] * conj(psiView[k]) * phiView[k];
    }
    
    // reduce and scatter such that energy is available on
    // all cores
    int count = 1; // send *one* integer
    
    Teuchos::Array<double_complex> sendBuff ( count );
    sendBuff[0] = localSum;
    
    Teuchos::Array<double_complex> recvBuff ( count );
    
    Teuchos::reduceAll ( * ( psi_.getMap()->getComm() ),
                         Teuchos::REDUCE_SUM,
                         count,
                         sendBuff.getRawPtr(),
                         recvBuff.getRawPtr()
                       );
                                 
    return recvBuff[0];
}
// ============================================================================
double
Ginla::FVM::State::
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
Ginla::FVM::State::
update( const double                  alpha,
        const Ginla::State::Virtual & b,
        const double                  beta
      )
{
  psi_.update( alpha, *(b.getPsi()), beta );
  chi_ = alpha*chi_ + beta * b.getChi();
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
Ginla::FVM::State::
getVorticity () const
{
//   TEST_FOR_EXCEPTION( true,
//                       std::logic_error,
//                       "Method \"getVorticity()\" not yet implemented."
//                     );

  return 0;
}
// ============================================================================