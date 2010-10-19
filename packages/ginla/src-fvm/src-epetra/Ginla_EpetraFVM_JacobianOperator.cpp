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

#include "Ginla_EpetraFVM_JacobianOperator.h"

#include <Teuchos_ArrayRCP.hpp>

// =============================================================================
Ginla::EpetraFVM::JacobianOperator::
JacobianOperator( const Teuchos::RCP<VIO::EpetraMesh::Mesh>                   & mesh,
                  const Teuchos::RCP<Ginla::EpetraFVM::KineticEnergyOperator> & keo
                ):
        useTranspose_( false ),
        comm_( mesh->getNodesMap()->Comm() ),
        mesh_( mesh ),
        keo_( keo ),
        currentX_ ( Teuchos::null ),
        temperature_( 0.0 )
{
}
// =============================================================================
Ginla::EpetraFVM::JacobianOperator::
~JacobianOperator ()
{
}
// =============================================================================
int
Ginla::EpetraFVM::JacobianOperator::
SetUseTranspose( bool UseTranspose )
{
    useTranspose_ = UseTranspose;
    return 0;
}
// =============================================================================
int
Ginla::EpetraFVM::JacobianOperator::
Apply ( const Epetra_MultiVector & X,
              Epetra_MultiVector & Y
      ) const
{
    // Add the terms corresponding to the nonlinear terms.
    // A = K + I * ( 1 - 2*|psi|^2 )
    // B = diag( -psi^2 )

    // K*psi
    TEUCHOS_ASSERT_EQUALITY( 0, keo_->Apply( X, Y ) )

//    TEUCHOS_ASSERT_EQUALITY( 0, 0 );

    for ( int vec=0; vec<X.NumVectors(); vec++ )
    {
        for ( int k=0; k<X.MyLength(); k++ )
        {
            // terms corresponding to  I * ( (1-T) - 2*|psi|^2 )
            double alpha = (1.0-temperature_)
                           - 2.0 * ( (*currentX_)[2*k]*(*currentX_)[2*k] + (*currentX_)[2*k+1]*(*currentX_)[2*k+1] );
            // real part
            Y.SumIntoMyValue( 2*k, vec, alpha * X[vec][2*k] );
            // imaginary part
            Y.SumIntoMyValue( 2*k+1, vec, alpha * X[vec][2*k+1] );

            // terms corresponding to  B = diag( -psi^2 )
            // Re(phi^2)
            double rePhiSquare = (*currentX_)[2*k]*(*currentX_)[2*k] - (*currentX_)[2*k+1]*(*currentX_)[2*k+1];
            // Im(phi^2)
            double imPhiSquare = 2.0*(*currentX_)[2*k]*(*currentX_)[2*k+1];
            // real part
            Y.SumIntoMyValue( 2*k,   vec, - rePhiSquare * X[vec][2*k] + imPhiSquare * X[vec][2*k+1] );
            // imaginary part
            Y.SumIntoMyValue( 2*k+1, vec, - imPhiSquare * X[vec][2*k] - rePhiSquare * X[vec][2*k+1] );
        }
    }

    return 0;
}
// =============================================================================
int
Ginla::EpetraFVM::JacobianOperator::
ApplyInverse ( const Epetra_MultiVector & X,
                     Epetra_MultiVector & Y
             ) const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Not yet implemented." );
    return -1;
}
// =============================================================================
double
Ginla::EpetraFVM::JacobianOperator::
NormInf () const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Not yet implemented." );
    return 0.0;
}
// =============================================================================
const char *
Ginla::EpetraFVM::JacobianOperator::
Label () const
{
    return "Jacobian operator for Ginzburg--Landau";
}
// =============================================================================
bool
Ginla::EpetraFVM::JacobianOperator::
UseTranspose () const
{
    return useTranspose_;
}
// =============================================================================
bool
Ginla::EpetraFVM::JacobianOperator::
HasNormInf () const
{
    return false;
}
// =============================================================================
const Epetra_Comm &
Ginla::EpetraFVM::JacobianOperator::
Comm () const
{
    return comm_;
}
// =============================================================================
const Epetra_Map &
Ginla::EpetraFVM::JacobianOperator::
OperatorDomainMap () const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Not yet implemented." );
}
// =============================================================================
const Epetra_Map &
Ginla::EpetraFVM::JacobianOperator::
OperatorRangeMap () const
{
    TEST_FOR_EXCEPTION( true,
                        std::logic_error,
                        "Not yet implemented." );
}
// =============================================================================
void
Ginla::EpetraFVM::JacobianOperator::
setParameters( const double mu,
               const Teuchos::Tuple<double,3> & scaling,
               const double temperature
             )
{
    keo_->setParameters( mu, scaling );
    temperature_ = temperature;
    return;
}
// =============================================================================
void
Ginla::EpetraFVM::JacobianOperator::
setCurrentX( const Teuchos::RCP<const Epetra_Vector> & currentX )
{
    currentX_ = currentX;
    return;
}
// =============================================================================
