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

#include "Ginla_FVM_EpetraModelEvaluator.h"

#include "Komplex2_LinearProblem.h"
#include "Komplex2_DoubleMatrix.h"

#include "Ginla_FVM_State.h"
#include "Ginla_State_Virtual.h"
#include "Ginla_MagneticVectorPotential_Virtual.h"

#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialDenseMatrix.h>

// ============================================================================
Ginla::FVM::ModelEvaluator::
ModelEvaluator ( const Teuchos::RCP<VIO::Mesh::Mesh>                         & mesh,
                 const Teuchos::ParameterList                                & problemParams,
                 const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp,
                 const Teuchos::RCP<Komplex2::LinearProblem>                 & komplex,
                 const Teuchos::RCP<Ginla::FVM::State>                       & initialState
               ) :
        mesh_ ( mesh ),
        komplex_ ( komplex ),
        mvp_( mvp ),
        x_(  Teuchos::null ),
        firstTime_ ( true ),
        numParams_( 1 ),
        p_map_( Teuchos::null ),
        p_init_( Teuchos::null ),
        p_names_( Teuchos::null ),
        p_current_( Teuchos::null )
//        jacobianOperator_( Teuchos::rcp( new Komplex2::DoubleMatrix( komplex_->getComplexMap(),
//                                                                     komplex_->getComplexMap()
//                                                                   )
//                                       )
//                         )
{
  this->setupParameters_( problemParams );

  // prepare initial guess
//   Teuchos::RCP<Ginla::FVM::State> initialState =
//       Teuchos::rcp( new Ginla::FVM::State( komplex_->getComplexMap(), mesh_) );

//   initialState->getPsiNonConst()->putScalar( double_complex(1.0,0.0) );
//   initialState = state;
//   initialState->getPsiNonConst()->randomize();

  x_ = this->createSystemVector( *initialState );

  return;
}
// ============================================================================
Ginla::FVM::ModelEvaluator::
~ModelEvaluator()
{
}
// ============================================================================
void
Ginla::FVM::ModelEvaluator::
setupParameters_( const Teuchos::ParameterList & params )
{
  p_names_ = Teuchos::rcp( new Teuchos::Array<std::string>() );
  p_names_->append( "mu" );
  p_names_->append( "scaling" );
  p_names_->append( "scaling x" );
  p_names_->append( "scaling y" );
  p_names_->append( "scaling z" );
  p_names_->append( "temperature" );

  // these are local variables
  Teuchos::RCP<Teuchos::Array<double> > p_default_values = Teuchos::rcp( new Teuchos::Array<double>() );
  p_default_values->append( 0.0 );
  p_default_values->append( 1.0 );
  p_default_values->append( 1.0 );
  p_default_values->append( 1.0 );
  p_default_values->append( 1.0 );
  p_default_values->append( 0.0 );

  // setup parameter map
  numParams_ = p_names_->length();
  const Epetra_Comm & comm = komplex_->getRealMap()->Comm();
  p_map_ = Teuchos::rcp( new Epetra_LocalMap( numParams_,
                                              0,
                                              comm
                                            )
                       );

  // run through the p_init and fill it with either what's given in params,
  // or the default value
  p_init_ = Teuchos::rcp( new Epetra_Vector(*p_map_) );
  for ( int k=0; k<numParams_; k++ )
      if ( params.isParameter( (*p_names_)[k] ) )
          (*p_init_)[k] = params.get<double>( (*p_names_)[k] );
      else
      {
          (*p_init_)[k] = (*p_default_values)[k];
          std::cerr << "Parameter \"" << (*p_names_)[k]
                    << "\" initialized with default value \""
                    << (*p_default_values)[k] << "\"."
                    << std::endl;
      }

  // TODO warn if there are unused entries in params
//  for ( Teuchos::ParameterList::ConstIterator k=params.begin(); k!=params.end(); ++k )
//  {
//      if ( !p_names_)
//      std::cout << "Parameter " << params.name(k) << std::endl;
//  }

  // also initialize p_current_
  p_current_ = Teuchos::rcp( new Epetra_Vector(*p_map_) );

  return;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::FVM::ModelEvaluator::
get_x_map() const
{
  return komplex_->getRealMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::FVM::ModelEvaluator::
get_f_map() const
{
  return komplex_->getRealMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::FVM::ModelEvaluator::
get_x_init () const
{
  TEUCHOS_ASSERT( !x_.is_null() );
  return x_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::FVM::ModelEvaluator::
get_p_init ( int l ) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_init_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::FVM::ModelEvaluator::
get_p_map(int l) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_map_;
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
Ginla::FVM::ModelEvaluator::
get_p_names (int l) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_names_;
}
// ============================================================================
Teuchos::RCP<Epetra_Operator>
Ginla::FVM::ModelEvaluator::
create_W() const
{
  TEUCHOS_ASSERT( !komplex_.is_null() );
  return komplex_->getMatrix();
}
// ============================================================================
EpetraExt::ModelEvaluator::InArgs
Ginla::FVM::ModelEvaluator::
createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  inArgs.setModelEvalDescription( "FVM Ginzburg-Landau" );

  // TODO is this actually correct?
  inArgs.set_Np( numParams_ );

  inArgs.setSupports( IN_ARG_x, true );

//   // for shifted matrix
//   inArgs.setSupports( IN_ARG_alpha, true );
//   inArgs.setSupports( IN_ARG_beta, true );

  return inArgs;
}
// ============================================================================
EpetraExt::ModelEvaluator::OutArgs
Ginla::FVM::ModelEvaluator::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  outArgs.setModelEvalDescription( "FVM Ginzburg-Landau" );

  outArgs.set_Np_Ng( 1, 0 ); // one parameter vector, no objective function

  // support derivatives with respect to all parameters
  outArgs.setSupports( OUT_ARG_DfDp,
                       0,
                       DerivativeSupport(DERIV_MV_BY_COL)
                     );

  outArgs.setSupports( OUT_ARG_f, true );
  outArgs.setSupports( OUT_ARG_W, true );

  outArgs.set_W_properties( DerivativeProperties( DERIV_LINEARITY_UNKNOWN, // DERIV_LINEARITY_NONCONST
                                                  DERIV_RANK_UNKNOWN, // DERIV_RANK_FULL, DERIV_RANK_DEFICIENT
                                                  false // supportsAdjoint
                                                )
                          );

  return outArgs;
}
// ============================================================================
void
Ginla::FVM::ModelEvaluator::
evalModel( const InArgs  & inArgs,
           const OutArgs & outArgs
         ) const
{
//   double alpha = inArgs.get_alpha();
//   double beta  = inArgs.get_beta();

  const Teuchos::RCP<const Epetra_Vector> & x_in = inArgs.get_x();

  // get input arguments and make sure they are all right
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  TEUCHOS_ASSERT( !p_in.is_null() );
  for ( int k=0; k<p_in->MyLength(); k++ )
      TEUCHOS_ASSERT( !isnan( (*p_in)[k] ) );

  const double mu = (*p_in)[0];
  const double scaling = (*p_in)[1];
  const Teuchos::Tuple<double,3> scalingX = Teuchos::tuple( (*p_in)[2],
                                                            (*p_in)[3],
                                                            (*p_in)[4]
                                                          );
  const double temperature = (*p_in)[5];

  // p_current_ is used in getParameters.
  // Setting p_current_=p_in here is really a somewhat arbitrary choice.
  // The rationale is that the p-values which were last
  // used here are the "current" parameter values,
  // but there's actually no guarantee for it:
  // The evaluator could have been used for *anything.
  // Anyway, current_p_ is only used in this->getParameters().
  *p_current_ = *p_in;

  Teuchos::Tuple<double,3> scalingCombined = Teuchos::tuple( scaling * scalingX[0],
                                                             scaling * scalingX[1],
                                                             scaling * scalingX[2]
                                                           );

  // compute F
  const Teuchos::RCP<Epetra_Vector> f_out = outArgs.get_f();
  if ( !f_out.is_null() )
      this->computeF_( *x_in, mu, scalingCombined, temperature, *f_out );

  // fill jacobian
  const Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();
  if( !W_out.is_null() )
  {
//      Teuchos::RCP<Epetra_CrsMatrix> W_out_crs =
//          Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>( W_out, true );
      this->computeJacobian_( *x_in, mu, scalingCombined, temperature, *W_out_crs );

//       W_out_crs->Scale( beta );
//
//       double diag = -alpha;
//       for ( int i=0; i<x_in->MyLength(); i++ )
//           W_out_crs->SumIntoMyValues( i, 1, &diag, &i );
  }

  return;
}
// ============================================================================
void
Ginla::FVM::ModelEvaluator::
computeF_ ( const Epetra_Vector            & x,
            const double                     mu,
            const Teuchos::Tuple<double,3> & scaling,
            const double                     temperature,
            Epetra_Vector                  & FVec
          ) const
{
  // compute FVec = K*x
  keo_->setParameters( mu, scaling );
  keo_->Apply( x, FVec );

  // add the nonlinear part (mass lumping)
  TEUCHOS_ASSERT( FVec.Map().SameAs( x.Map() ) );

  for ( int k=0; k<controlVolumes.MyLength(); k++ )
  {
      // Do the equivalent of
      //   res[k] -= controlVolumes[k] * psi[k] * ( (1.0-temperature) - std::norm(psi[k]) );
      double alpha = - controlVolumes[k]
                     * ( (1.0-temperature) - x[2*k]*x[2*k] - x[2*k+1]*x[2*k+1] );
      // real part
      FVec.SumIntoMyValue( 2*k,
                           0,
                           alpha * x[2*k] );
      // imaginary part
      FVec.SumIntoMyValue( 2*k+1,
                           0,
                           alpha * x[2*k+1] );
  }

  return;
}
// ============================================================================
void
Ginla::FVM::ModelEvaluator::
computeJacobian_ ( const Epetra_Vector            & x,
                   const double                     mu,
                   const Teuchos::Tuple<double,3> & scaling,
                   const double                     temperature,
                   Epetra_CrsMatrix               & Jac
                 ) const
{
  jacobianOperator_->setParameters( mu, scaling, temperature );
  jacobianOperator_->setCurrentX( x );

  return;
}
// ============================================================================
Teuchos::RCP<Ginla::State::Virtual>
Ginla::FVM::ModelEvaluator::
createState( const Epetra_Vector & x ) const
{
    const Teuchos::RCP<ComplexVector> psi = komplex_->real2complex ( x );
    return Teuchos::rcp( new Ginla::FVM::State( psi, mesh_ ) );
}
// =============================================================================
Teuchos::RCP<Epetra_Vector>
Ginla::FVM::ModelEvaluator::
createSystemVector(  const Ginla::State::Virtual & state ) const
{
    return komplex_->complex2real ( state.getPsi() );
}
// =============================================================================
void
Ginla::FVM::ModelEvaluator::
createSystemVector( const Ginla::State::Virtual & state,
                          Epetra_Vector         & x
                  ) const
{
  komplex_->complex2real( *state.getPsi(), x );
  return;
}
// =============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::FVM::ModelEvaluator::
getParameters() const
{
  // construct a LOCA::ParameterVector of the parameters
  Teuchos::RCP<LOCA::ParameterVector> p =
      Teuchos::rcp( new LOCA::ParameterVector() );

  TEUCHOS_ASSERT( !p_names_.is_null() );
  TEUCHOS_ASSERT( !p_current_.is_null() );
  for ( int k=0; k<numParams_; k++ )
      p->addParameter( (*p_names_)[k], (*p_current_)[k] );

  return p;
}
// =============================================================================
