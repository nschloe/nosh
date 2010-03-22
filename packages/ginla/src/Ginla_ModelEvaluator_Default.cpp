/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010 Nico Sch\"omer

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

#include "Ginla_ModelEvaluator_Default.h"

#include <Epetra_Map.h>

// ============================================================================
Ginla::ModelEvaluator::Default::
Default ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
          const Teuchos::RCP<Ginla::Komplex>           & komplex
        ) :
        glOperator_ ( glOperator ),
        komplex_ ( komplex ),
        x_(  Teuchos::rcp( new Epetra_Vector( *komplex_->getRealMap(), true ) ) ),
        firstTime_ ( true )
{
//   x_->Random();
  x_->PutScalar( 1.0 );
}
// ============================================================================
Ginla::ModelEvaluator::Default::
Default ( const Teuchos::RCP<Ginla::Operator::Virtual> & glOperator,
          const Teuchos::RCP<Ginla::Komplex>           & komplex,
          const Teuchos::RCP<const ComplexVector>      & psi
        ) :
        glOperator_ ( glOperator ),
        komplex_ ( komplex ),
        x_( komplex_->complex2real( psi ) ),
        firstTime_ ( true )
{
  // make sure the maps are compatible
  TEUCHOS_ASSERT( psi->getMap()->isSameAs( *komplex->getComplexMap()) );
}
// ============================================================================
Ginla::ModelEvaluator::Default::
~Default()
{
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::ModelEvaluator::Default::
get_x_map() const
{
  return komplex_->getRealMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::ModelEvaluator::Default::
get_f_map() const
{
  return komplex_->getRealMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::ModelEvaluator::Default::
get_x_init () const
{
  return x_;
}
// ============================================================================
Teuchos::RCP<Epetra_Operator>
Ginla::ModelEvaluator::Default::
create_W() const
{
  return komplex_->getMatrix();
}
// ============================================================================
EpetraExt::ModelEvaluator::InArgs
Ginla::ModelEvaluator::Default::
createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription( "Ginzburg--Landau extreme type-II on a square" );
  inArgs.setSupports( IN_ARG_x,
                      true // supports
                    );
  return inArgs;
}
// ============================================================================
EpetraExt::ModelEvaluator::OutArgs
Ginla::ModelEvaluator::Default::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription("Ginzburg--Landau extreme type-II on a square");
  outArgs.setSupports( OUT_ARG_f, true );
  outArgs.setSupports( OUT_ARG_W, true );
  outArgs.set_W_properties( DerivativeProperties( DERIV_LINEARITY_NONCONST,
                                                  DERIV_RANK_FULL,
                                                  false // supportsAdjoint
                                                )
                          );
  return outArgs;
}
// ============================================================================
void
Ginla::ModelEvaluator::Default::
evalModel( const InArgs  & inArgs, 
           const OutArgs & outArgs
         ) const
{
//   std::cout << "\n\n>>>>>>>>>>>>>>>>>>>>> evalModel\n";
//   std::cout << outArgs.get_W().is_null() << " " << outArgs.get_f().is_null()
//             << std::endl;
  
  const Teuchos::RCP<const Epetra_Vector> & x = inArgs.get_x();
  
  Teuchos::RCP<Epetra_Vector>   f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();
  
//   std::cout << "norm at 0: ";
//   double alpha[1];
//   x->Norm2( alpha );
//   std::cout << alpha[0] << std::endl;
  
  // compute F
  if (!f_out.is_null())
     this->computeF_( *x, *f_out ); 
  
//   std::cout << "after: ";
//   outArgs.get_f()->Norm2( alpha );
//   std::cout << alpha[0] << std::endl;

  if( !W_out.is_null() )
  {
      // fill jacobian
      Teuchos::RCP<Epetra_RowMatrix> tmp =
             Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(W_out);
      this->computeJacobian_( *x, *W_out );
  }

  return;
}
// ============================================================================
void
Ginla::ModelEvaluator::Default::
computeF_ ( const Epetra_Vector & x,
            Epetra_Vector       & FVec
          ) const
{
    // ------------------------------------------------------------------------
    // convert from x to psi
    const Teuchos::RCP<ComplexVector> psi = komplex_->real2complex ( x );
    // ------------------------------------------------------------------------
    // compute the GL residual
    // setup output vector with the same map as psi
    Teuchos::RCP<ComplexVector> res =
        Teuchos::rcp ( new ComplexVector ( psi->getMap(), true ) );

    // TODO not really necessary?
    glOperator_->updatePsi ( psi );

    Teuchos::ArrayRCP<double_complex> resView = res->get1dViewNonConst();
    // loop over the nodes
    for ( unsigned int k=0; k<psi->getLocalLength(); k++ )
    {
        int globalIndex = psi->getMap()->getGlobalElement ( k );
        resView[k] = glOperator_->getEntry ( globalIndex );

//         if ( !perturbation_.is_null() )
//             resView[k] += perturbation_->computePerturbation ( globalIndex );
    }
    // ------------------------------------------------------------------------
    // TODO Avoid this explicit copy?
    // transform back to fully real equation
    FVec = * ( komplex_->complex2real ( *res ) );
    // ------------------------------------------------------------------------
    
    return;
}
// ============================================================================
void
Ginla::ModelEvaluator::Default::
computeJacobian_ ( const Epetra_Vector & x,
                   Epetra_Operator     & Jac
                 ) const
{
    Teuchos::Array<int> indicesA, indicesB;
    Teuchos::Array<double_complex> valuesA, valuesB;

    komplex_->zeroOutMatrix();

    Teuchos::RCP<ComplexVector> psi = komplex_->real2complex ( x );

    // update to the latest psi vector before retrieving the Jacobian
    glOperator_->updatePsi ( psi );
    
    // loop over the rows and fill the matrix
    int numMyElements = komplex_->getComplexMap()->getNodeNumElements();
    for ( int row = 0; row < numMyElements; row++ )
    {
        int globalRow = komplex_->getComplexMap()->getGlobalElement ( row );
        // get the values from the operator
        glOperator_->getJacobianRow ( globalRow,
                                      indicesA, valuesA,
                                      indicesB, valuesB );
        // ... and fill them into glKomplex_
        komplex_->updateRow ( globalRow,
                              indicesA, valuesA,
                              indicesB, valuesB,
                              firstTime_ );
    }

    if ( firstTime_ )
    {
        komplex_->finalizeMatrix();
        firstTime_ = false;
    }

    return;
}
// ============================================================================
Teuchos::RCP<const Ginla::Komplex>
Ginla::ModelEvaluator::Default::
getKomplex() const
{
    return komplex_; 
}
// ============================================================================