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
#include "Ginla_MagneticVectorPotential_Virtual.h"
#include "Ginla_EpetraFVM_KeoFactory.h"
#include "Ginla_EpetraFVM_KeoPreconditioner.h"
#include "Ginla_EpetraFVM_StkMesh.h"

#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialDenseMatrix.h>

// ============================================================================
Ginla::EpetraFVM::ModelEvaluator::
ModelEvaluator ( const Teuchos::RCP<Ginla::EpetraFVM::StkMesh>               & mesh,
                 const Teuchos::ParameterList                                & problemParams,
                 const Teuchos::RCP<const Epetra_Vector>                     & thickness,
                 const Teuchos::RCP<Ginla::MagneticVectorPotential::Virtual> & mvp,
                 const Teuchos::RCP<Ginla::EpetraFVM::State>                 & initialState
               ) :
        mesh_ ( mesh ),
        thickness_( thickness ),
        x_( initialState->getPsiNonConst() ),
        numParams_( 1 ),
        p_map_( Teuchos::null ),
        p_init_( Teuchos::null ),
        p_names_( Teuchos::null ),
        p_current_( Teuchos::null ),
        mvp_( mvp ),
        keoFactory_( Teuchos::rcp( new Ginla::EpetraFVM::KeoFactory( mesh, thickness, mvp ) ) )
{
  this->setupParameters_( problemParams );

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
  p_names_->append( "phi" );
  p_names_->append( "theta" );
  p_names_->append( "scaling" );
  p_names_->append( "scaling x" );
  p_names_->append( "scaling y" );
  p_names_->append( "scaling z" );
  p_names_->append( "temperature" );

  // these are local variables
  Teuchos::RCP<Teuchos::Array<double> > p_default_values = Teuchos::rcp( new Teuchos::Array<double>() );
  p_default_values->append( 0.0 ); // mu
  p_default_values->append( 0.0 ); // phi
  p_default_values->append( 0.0 ); // theta
  p_default_values->append( 1.0 ); // scaling
  p_default_values->append( 1.0 ); // scaling x
  p_default_values->append( 1.0 ); // scaling y
  p_default_values->append( 1.0 ); // scaling z
  p_default_values->append( 0.0 ); // temperature

  // setup parameter map
  numParams_ = p_names_->length();
  TEUCHOS_ASSERT( !mesh_.is_null() );
  const Epetra_Comm & comm = mesh_->getComm();
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
          if ( comm.MyPID() == 0 )
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
  return mesh_->getComplexMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Ginla::EpetraFVM::ModelEvaluator::
get_f_map() const
{
    TEUCHOS_ASSERT( !mesh_.is_null() );
    return mesh_->getComplexMap();
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
  return Teuchos::rcp( new Ginla::EpetraFVM::JacobianOperator( mesh_, thickness_, mvp_ ) );
}
// =============================================================================
Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Ginla::EpetraFVM::ModelEvaluator::
create_WPrec() const
{
  Teuchos::RCP<Epetra_Operator> keoPrec =
          Teuchos::rcp( new Ginla::EpetraFVM::KeoPreconditioner( mesh_, thickness_, mvp_ ) );
  // bool is answer to: "Prec is already inverted?"
  return Teuchos::rcp( new EpetraExt::ModelEvaluator::Preconditioner( keoPrec, false ) );
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

//   double norm1;
//   x_in->Norm1( &norm1 );
//   std::cout << "\n\n\t\tnorm1 = " << norm1 << "\n\n" << std::endl;

  // get input arguments and make sure they are all right
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  TEUCHOS_ASSERT( !p_in.is_null() );
  for ( int k=0; k<p_in->MyLength(); k++ )
      TEUCHOS_ASSERT( !isnan( (*p_in)[k] ) );

  // p_current_ is used in getParameters.
  // Setting p_current_=p_in here is really a somewhat arbitrary choice.
  // The rationale is that the p-values which were last
  // used here are the "current" parameter values,
  // but there's actually no guarantee for it:
  // The evaluator could have been used for *anything.
  // Anyway, current_p_ is only used in this->getParameters().
  *p_current_ = *p_in;

  Teuchos::RCP<LOCA::ParameterVector> mvpParams = this->getParameters();

  const double scaling = (*p_in)[3];
  const Teuchos::Tuple<double,3> scalingX = Teuchos::tuple( (*p_in)[4],
                                                            (*p_in)[5],
                                                            (*p_in)[6]
                                                          );
  const double temperature = (*p_in)[7];

  Teuchos::Tuple<double,3> scalingCombined = Teuchos::tuple( scaling * scalingX[0],
                                                             scaling * scalingX[1],
                                                             scaling * scalingX[2]
                                                           );

  // compute F
  const Teuchos::RCP<Epetra_Vector> f_out = outArgs.get_f();
  if ( !f_out.is_null() )
  {
      this->computeF_( *x_in, mvpParams, scalingCombined, temperature, *f_out );
  }

  // fill jacobian
  const Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();
  if( !W_out.is_null() )
  {
      Teuchos::RCP<Ginla::EpetraFVM::JacobianOperator> jac =
          Teuchos::rcp_dynamic_cast<Ginla::EpetraFVM::JacobianOperator>( W_out, true );
      jac->rebuild( mvpParams,
                    scalingCombined,
                    temperature,
                    x_in
                  );
  }

  // fill preconditioner
  const Teuchos::RCP<Epetra_Operator> WPrec_out = outArgs.get_WPrec();
  if( !WPrec_out.is_null() )
  {
      Teuchos::RCP<Ginla::EpetraFVM::KeoPreconditioner> keoPrec =
          Teuchos::rcp_dynamic_cast<Ginla::EpetraFVM::KeoPreconditioner>( WPrec_out, true );
      keoPrec->rebuild( mvpParams,
                        scalingCombined );
  }

  return;
}
// ============================================================================
void
Ginla::EpetraFVM::ModelEvaluator::
computeF_ ( const Epetra_Vector                             & x,
            const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams,
            const Teuchos::Tuple<double,3>                  & scaling,
            const double                                      temperature,
            Epetra_Vector                                   & FVec
          ) const
{
  // build the KEO
  Epetra_FECrsMatrix keoMatrix( Copy, keoFactory_->buildKeoGraph() );
  keoFactory_->updateParameters( mvpParams, scaling );
  keoFactory_->buildKeo( keoMatrix );

  // compute FVec = K*x
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix.Apply( x, FVec ) );

  // add the nonlinear part (mass lumping)
  TEUCHOS_ASSERT( FVec.Map().SameAs( x.Map() ) );

  const Epetra_Vector & controlVolumes = *(mesh_->getControlVolumes());

  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY( 2*controlVolumes.MyLength(), x.MyLength() );

  for ( int k=0; k<controlVolumes.MyLength(); k++ )
  {
      // In principle, mass lumping here suggests to take
      //
      //   \int_{control volume}  thickness * f(psi).
      //
      // with f(psi) = psi[k] * ( (1.0-temperature) - std::norm(psi[k]) ).
      // (a) A possible approximation for this is
      //
      //        |control volume| * average(thicknesses) * f(psi(x_k)).
      //
      //     This is loosely derived from the midpoint quadrature rule for
      //     triangles, i.e.,
      //
      //      \int_{triangle} f(x) ~= |triangle| *  \sum_{edge midpoint} 1/3 * f(midpoint).
      //
      //     so all the values at the midpoints have the same weight (independent of
      //     whether the edge is long or short). This is then "generalized to
      //
      //      \int_{triangle} f(x)*a(x) ~= |triangle| *  \sum_{edge midpoint} 1/3 * f(midpoint)*a(midpoint),
      //
      //     or, as f(midpoint) is not available,
      //
      //      \int_{triangle} f(x)*a(x) ~= |triangle| * f(center of gravity)  \sum_{edge midpoint} 1/3 * a(midpoint).
      //
      //     For general polynomals, this is then the above expression.
      //     Hence, do the equivalent of
      //
      //       res[k] += controlVolumes[k] * average(thicknesses) * psi[k] * ( (1.0-temperature) - std::norm(psi[k]) );
      //
      // (b) Another possible approximation is 
      //
      //        |control volume| * thickness(x_k) * f(psi(x_k))
      //
      //     as suggested by mass lumping. This works if if thickness(x_k) is available.
      //
      double alpha = controlVolumes[k] * (*thickness_)[k]
                     * ( (1.0-temperature) - x[2*k]*x[2*k] - x[2*k+1]*x[2*k+1] );

      // real part
      TEUCHOS_ASSERT_EQUALITY( 0, FVec.SumIntoMyValue( 2*k,
                                                       0,
                                                       alpha * x[2*k] )
                             );
      // imaginary part
      TEUCHOS_ASSERT_EQUALITY( 0, FVec.SumIntoMyValue( 2*k+1,
                                                       0,
                                                       alpha * x[2*k+1] )
                             );
  }

  return;
}
// ============================================================================
Teuchos::RCP<Ginla::EpetraFVM::State>
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
