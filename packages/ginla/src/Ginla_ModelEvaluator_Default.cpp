/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2010  Nico Sch\"omer

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

#include <Epetra_LocalMap.h>
#include <Epetra_CrsMatrix.h>

#include "Ginla_Operator_Virtual.h"

// ============================================================================
Ginla::ModelEvaluator::Default::
Default ( const Teuchos::RCP<Ginla::Operator::Virtual>      & glOperator,
          const Teuchos::RCP<Ginla::Komplex::LinearProblem> & komplex
        ) :
        glOperator_ ( glOperator ),
        komplex_ ( komplex ),
        x_(  Teuchos::rcp( new Epetra_Vector( *komplex_->getRealMap(), true ) ) ),
        firstTime_ ( true ),
        numParameters_( 1 )
{ 
//   x_->Random();

  Teuchos::RCP<ComplexVector> psi =
      Teuchos::rcp( new ComplexVector( komplex_->getComplexMap() ) );
  psi->putScalar( double_complex(1.0,0.0) );
  *x_ = *(komplex_->complex2real( psi ));
  
  // set up parameters
  p_map_ = Teuchos::rcp(new Epetra_LocalMap( numParameters_,
                                             0,
                                             komplex_->getRealMap()->Comm() ) );

  p_init_ = Teuchos::rcp(new Epetra_Vector(*p_map_));
  for (int i=0; i<numParameters_; i++)
      (*p_init_)[i]= 0.0;
  
  Teuchos::Tuple<std::string,1> t = Teuchos::tuple<std::string>( "H0" );
  p_names_ = Teuchos::rcp( new Teuchos::Array<std::string>( t ) );

  return;
}
// ============================================================================
Ginla::ModelEvaluator::Default::
Default ( const Teuchos::RCP<Ginla::Operator::Virtual>      & glOperator,
          const Teuchos::RCP<Ginla::Komplex::LinearProblem> & komplex,
          const Teuchos::RCP<const Ginla::State>            & state
        ) :
        glOperator_ ( glOperator ),
        komplex_ ( komplex ),
        x_( this->createSystemVector_( *state ) ),
        firstTime_ ( true ),
        numParameters_( 1 )
{
  // set up parameters
  p_map_ = Teuchos::rcp(new Epetra_LocalMap( numParameters_,
                                             0,
                                             komplex_->getRealMap()->Comm() ) );

  p_init_ = Teuchos::rcp(new Epetra_Vector(*p_map_));
  for (int i=0; i<numParameters_; i++)
      (*p_init_)[i]= 0.0;
  
  // make sure the maps are compatible
  TEUCHOS_ASSERT( state->getPsi()->getMap()->isSameAs( *komplex->getComplexMap()) );
  
  Teuchos::Tuple<std::string,1> t = Teuchos::tuple<std::string>( "H0" );
  p_names_ = Teuchos::rcp( new Teuchos::Array<std::string>( t ) );
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
Teuchos::RCP<const Epetra_Vector>
Ginla::ModelEvaluator::Default::
get_p_init ( int l ) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_init_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::ModelEvaluator::Default::
get_p_map(int l) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_map_;
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
Ginla::ModelEvaluator::Default::
get_p_names (int l) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_names_;
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
  
  inArgs.set_Np( 1 );
  
  inArgs.setSupports( IN_ARG_x, true );
  
  return inArgs;
}
// ============================================================================
EpetraExt::ModelEvaluator::OutArgs
Ginla::ModelEvaluator::Default::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  
  outArgs.setModelEvalDescription( "Ginzburg--Landau extreme type-II on a square" );
  
  outArgs.set_Np_Ng( 1 , 0 ); // return parameters p and solution g
  
  outArgs.setSupports( OUT_ARG_f, true );
  outArgs.setSupports( OUT_ARG_DfDp, 0, DerivativeSupport(DERIV_MV_BY_COL) );
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
  const Teuchos::RCP<const Epetra_Vector> & x_in = inArgs.get_x();
  
  Teuchos::RCP<Epetra_Vector>   f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();

  Teuchos::RCP<Epetra_MultiVector> dfdp_out;
  if ( outArgs.Np() > 0 )
      dfdp_out = outArgs.get_DfDp(0).getMultiVector();
  
  // Parse InArgs
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  
  TEUCHOS_ASSERT( !p_in.is_null() );
  
  // build a LOCA::ParameterVector list to feed glOperator_.
  LOCA::ParameterVector locaParams;
  for ( int k=0; k<p_names_->length(); k++ )
      locaParams.addParameter( (*p_names_)[k], (*p_in)[k] );
  glOperator_->setParameters( locaParams );
  
  // compute F
  if ( !f_out.is_null() )
      this->computeF_( *x_in, *f_out ); 

  // fill jacobian
  if( !W_out.is_null() )
      this->computeJacobian_( *x_in, *W_out );
  
  // dF / dH_0
  if ( !dfdp_out.is_null() )
      this->computeDFDh0_( *x_in, *((*dfdp_out)(0)) );

  return;
}
// ============================================================================
void
Ginla::ModelEvaluator::Default::
computeF_ ( const Epetra_Vector & x,
            Epetra_Vector       & FVec
          ) const
{
    // convert from x to state
    const Teuchos::RCP<const Ginla::State> state = this->createState( x );

    // compute the GL residual
    const Teuchos::RCP<const Ginla::State> res = glOperator_->getF( state );

    // TODO Avoid this explicit copy?
    // transform back to fully real equation
    FVec = *(this->createSystemVector_( *res ));

    return;
}
// ============================================================================
void
Ginla::ModelEvaluator::Default::
computeDFDh0_ ( const Epetra_Vector & x,
                Epetra_Vector       & FVec
              ) const
{
    // convert from x to state
    const Teuchos::RCP<const Ginla::State> state = this->createState( x );

    // compute the GL residual
    const Teuchos::RCP<const Ginla::State> res = glOperator_->getDFDh0( state );

    // TODO Avoid this explicit copy?
    // transform back to fully real equation
    FVec = *(this->createSystemVector_( *res ));

    return;
}
// ============================================================================
void
Ginla::ModelEvaluator::Default::
computeJacobian_ ( const Epetra_Vector & x,
                   Epetra_Operator     & Jac
                 ) const
{
    // In Ginla::Komplex, the values are merely sumIntoLocalValue'd,
    // so make sure we set this to zero from the start.
    komplex_->zeroOutMatrix();

    const Teuchos::RCP<const Ginla::State> state = this->createState( x );
    
    // create a real-valued matrix of the AB-Jacobian
    komplex_->update( glOperator_->getJacobian( state ), firstTime_ );

    if ( firstTime_ )
    {
        komplex_->finalizeMatrix();
        firstTime_ = false;
    }

    return;
}
// ============================================================================
Teuchos::RCP<const Ginla::Komplex::LinearProblem>
Ginla::ModelEvaluator::Default::
getKomplex() const
{
    return komplex_; 
}
// ============================================================================
Teuchos::RCP<Ginla::State>
Ginla::ModelEvaluator::Default::
createState( const Epetra_Vector & x ) const
{
    const Teuchos::RCP<ComplexVector> psi = komplex_->real2complex ( x );
    return Teuchos::rcp( new Ginla::State( psi, glOperator_->getGrid() ) );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::ModelEvaluator::Default::
createSystemVector_(  const Ginla::State & state ) const
{
    return komplex_->complex2real ( state.getPsi() );
}
// =============================================================================