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

#include "Ginla_EpetraFVM_ModelEvaluator.h"

#include "Ginla_EpetraFVM_State.h"
#include "Ginla_State_Virtual.h"
#include "Ginla_MagneticVectorPotential_Virtual.h"

#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialDenseMatrix.h>

// ============================================================================
Ginla::EpetraFVM::ModelEvaluator::
ModelEvaluator ( const Teuchos::RCP<VIO::EpetraMesh::Mesh>                   & mesh,
                 const Teuchos::ParameterList                                & problemParams,
                 const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp,
                 const Teuchos::RCP<Ginla::EpetraFVM::State>                 & initialState
               ) :
        mesh_ ( mesh ),
        x_( Teuchos::null ),
        numParams_( 1 ),
        p_map_( Teuchos::null ),
        p_init_( Teuchos::null ),
        p_names_( Teuchos::null ),
        p_current_( Teuchos::null ),
        keo_( Teuchos::rcp( new Ginla::EpetraFVM::KineticEnergyOperator( mesh, mvp ) ) )
{
  this->setupParameters_( problemParams );

  // prepare initial guess
  x_ = initialState->getPsiNonConst();

  return;
}
// ============================================================================
Ginla::EpetraFVM::ModelEvaluator::
~ModelEvaluator()
{
}
// ============================================================================
void
Ginla::EpetraFVM::ModelEvaluator::
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
  TEUCHOS_ASSERT( !mesh_.is_null() );
  const Epetra_Comm & comm = mesh_->getNodesMap()->Comm();
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
  p_current_ = Teuchos::rcp( new Epetra_Vector(*p_init_) );

  return;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::EpetraFVM::ModelEvaluator::
get_x_map() const
{
  TEUCHOS_ASSERT( !mesh_.is_null() );
  return mesh_->getComplexValuesMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::EpetraFVM::ModelEvaluator::
get_f_map() const
{
    TEUCHOS_ASSERT( !mesh_.is_null() );
    return mesh_->getComplexValuesMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::EpetraFVM::ModelEvaluator::
get_x_init () const
{
  TEUCHOS_ASSERT( !x_.is_null() );
  return x_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Ginla::EpetraFVM::ModelEvaluator::
get_p_init ( int l ) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_init_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::EpetraFVM::ModelEvaluator::
get_p_map( int l ) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_map_;
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
Ginla::EpetraFVM::ModelEvaluator::
get_p_names( int l ) const
{
  TEUCHOS_ASSERT_EQUALITY( 0, l );
  return p_names_;
}
// =============================================================================
Teuchos::RCP<Epetra_Operator>
Ginla::EpetraFVM::ModelEvaluator::
create_W() const
{
  return Teuchos::rcp( new Ginla::EpetraFVM::JacobianOperator( mesh_, keo_ ) );
}
// =============================================================================
Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Ginla::EpetraFVM::ModelEvaluator::
create_WPrec() const
{
  TEUCHOS_ASSERT( !keo_.is_null() );

  // bool is answer to: "Prec is already inverted?"
  return Teuchos::rcp( new EpetraExt::ModelEvaluator::Preconditioner( keo_,
                                                                      false
                                                                    )
                     );
}
// ============================================================================
EpetraExt::ModelEvaluator::InArgs
Ginla::EpetraFVM::ModelEvaluator::
createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  inArgs.setModelEvalDescription( "FVM Ginzburg-Landau" );

  // TODO is this actually correct?
  inArgs.set_Np( numParams_ );

  inArgs.setSupports( IN_ARG_x, true );

  // for shifted matrix
  inArgs.setSupports( IN_ARG_alpha, true );
  inArgs.setSupports( IN_ARG_beta, true );

  return inArgs;
}
// ============================================================================
EpetraExt::ModelEvaluator::OutArgs
Ginla::EpetraFVM::ModelEvaluator::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  outArgs.setModelEvalDescription( "FVM Ginzburg-Landau" );

  outArgs.set_Np_Ng( 1, 0 ); // one parameter vector, no objective function

  // support derivatives with respect to all parameters;
  // this is then handles by LOCA's finite differencing
  outArgs.setSupports( OUT_ARG_DfDp,
                       0,
                       DerivativeSupport(DERIV_MV_BY_COL)
                     );

  outArgs.setSupports( OUT_ARG_f, true );
  outArgs.setSupports( OUT_ARG_W, true );
  outArgs.set_W_properties( DerivativeProperties( DERIV_LINEARITY_UNKNOWN, // DERIV_LINEARITY_NONCONST
                                                  DERIV_RANK_DEFICIENT, // DERIV_RANK_FULL, DERIV_RANK_DEFICIENT
                                                  false // supportsAdjoint
                                                )
                          );

  outArgs.setSupports( OUT_ARG_WPrec, true );
  outArgs.set_WPrec_properties( DerivativeProperties( DERIV_LINEARITY_UNKNOWN,
                                                      DERIV_RANK_FULL,
                                                      false
                                                    )
                              );

  return outArgs;
}
// ============================================================================
void
Ginla::EpetraFVM::ModelEvaluator::
evalModel( const InArgs  & inArgs,
           const OutArgs & outArgs
         ) const
{
  const double alpha = inArgs.get_alpha();
  const double beta  = inArgs.get_beta();

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
      std::cout << "alpha, beta " << alpha << " " << beta << std::endl;

      Teuchos::RCP<Ginla::EpetraFVM::JacobianOperator> jac =
          Teuchos::rcp_dynamic_cast<Ginla::EpetraFVM::JacobianOperator>( W_out, true );

      jac->setParameters( mu, scalingCombined, temperature );
      jac->setShiftParameters( alpha, beta );
      jac->setCurrentX( x_in );
  }

  // fill preconditioner
  const Teuchos::RCP<Epetra_Operator> WPrec_out = outArgs.get_WPrec();
  if( !WPrec_out.is_null() )
  {
      Teuchos::RCP<Ginla::EpetraFVM::KineticEnergyOperator> keo =
          Teuchos::rcp_dynamic_cast<Ginla::EpetraFVM::KineticEnergyOperator>( WPrec_out, true );
      keo->setParameters( mu, scalingCombined );
  }

  return;
}
// ============================================================================
void
Ginla::EpetraFVM::ModelEvaluator::
computeF_ ( const Epetra_Vector            & x,
            const double                     mu,
            const Teuchos::Tuple<double,3> & scaling,
            const double                     temperature,
            Epetra_Vector                  & FVec
          ) const
{
  // compute FVec = K*x
  keo_->setParameters( mu, scaling );
  TEUCHOS_ASSERT_EQUALITY( 0, keo_->Apply( x, FVec ) );

  // add the nonlinear part (mass lumping)
  TEUCHOS_ASSERT( FVec.Map().SameAs( x.Map() ) );

  Epetra_Vector & controlVolumes = *(mesh_->getControlVolumes());

  for ( int k=0; k<controlVolumes.MyLength(); k++ )
  {
      // Do the equivalent of
      //   res[k] -= controlVolumes[k] * psi[k] * ( (1.0-temperature) - std::norm(psi[k]) );
      double alpha = controlVolumes[k]
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
Teuchos::RCP<Ginla::State::Virtual>
Ginla::EpetraFVM::ModelEvaluator::
createSavable( const Epetra_Vector & x ) const
{
    return Teuchos::rcp( new Ginla::EpetraFVM::State( x, mesh_ ) );
}
// =============================================================================
Teuchos::RCP<LOCA::ParameterVector>
Ginla::EpetraFVM::ModelEvaluator::
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
