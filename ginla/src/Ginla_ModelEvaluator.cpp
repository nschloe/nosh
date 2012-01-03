// @HEADER
//
//    <one line to give the program's name and a brief idea of what it does.>
//    Copyright (C) 2010, 2011  Nico Schl\"omer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// @HEADER

#include "Ginla_ModelEvaluator.hpp"

#include "Ginla_State.hpp"
#include "Ginla_MagneticVectorPotential.hpp"
#include "Ginla_KeoFactory.hpp"
#include "Ginla_KeoRegularized.hpp"
#include "Ginla_StkMesh.hpp"

#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialDenseMatrix.h>

#ifdef GINLA_TEUCHOS_TIME_MONITOR
  #include <Teuchos_TimeMonitor.hpp>
#endif

#include <Teuchos_VerboseObject.hpp>

namespace Ginla {
// ============================================================================
ModelEvaluator::
ModelEvaluator ( const Teuchos::RCP<Ginla::StkMesh>                 & mesh,
                 const Teuchos::ParameterList                       & problemParams,
                 const Teuchos::RCP<const Epetra_Vector>            & thickness,
                 const Teuchos::RCP<Ginla::MagneticVectorPotential> & mvp,
                 const Teuchos::RCP<Epetra_Vector>                  & initialX
               ) :
        mesh_ ( mesh ),
        thickness_( thickness ),
        x_( initialX ),
        numParams_( 0 ),
        p_map_( Teuchos::null ),
        p_init_( Teuchos::null ),
        p_names_( Teuchos::null ),
        p_current_( Teuchos::null ),
        mvp_( mvp ),
        keoFactory_( Teuchos::rcp( new Ginla::KeoFactory( mesh, thickness, mvp ) ) ),
#ifdef GINLA_TEUCHOS_TIME_MONITOR
        evalModelTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: ModelEvaluator::evalModel") ),
        computeFTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: ModelEvaluator::evalModel:compute F") ),
        computedFdpTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: ModelEvaluator::evalModel:compute dF/dp") ),
        fillJacobianTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: ModelEvaluator::evalModel:fill Jacobian") ),
        fillPreconditionerTime_( Teuchos::TimeMonitor::getNewTimer("Ginla: ModelEvaluator::fill preconditioner") ),
#endif
        out_( Teuchos::VerboseObjectBase::getDefaultOStream() )
{
  this->setupParameters_( problemParams );
  return;
}
// ============================================================================
ModelEvaluator::
~ModelEvaluator()
{
}
// ============================================================================
void
ModelEvaluator::
setupParameters_( const Teuchos::ParameterList & params )
{
  p_names_ = Teuchos::rcp( new Teuchos::Array<std::string>() );
  p_names_->append( "mu" );
  p_names_->append( "theta" );
  p_names_->append( "temperature" );

  numParams_ = p_names_->length();

  // these are local variables
  Teuchos::Array<double> p_default_values = Teuchos::Array<double>(numParams_);
  p_default_values[0] = 0.0; // mu
  p_default_values[1] = 0.0; // theta
  p_default_values[2] = 0.0; // temperature

  // setup parameter map
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
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
          (*p_init_)[k] = p_default_values[k];
          *out_ << "Parameter \"" << (*p_names_)[k]
                << "\" initialized with default value \""
                << p_default_values[k] << "\"."
                << std::endl;
      }

  // TODO warn if there are unused entries in params
//  for ( Teuchos::ParameterList::ConstIterator k=params.begin(); k!=params.end(); ++k )
//  {
//      if ( !p_names_)
//        *out_ << "Parameter " << params.name(k) << std::endl;
//  }

  // also initialize p_current_
  p_current_ = Teuchos::rcp( new Epetra_Vector(*p_init_) );

  return;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
ModelEvaluator::
get_x_map() const
{
  // It is a bit of an assumption that x_ actually has this map, but
  // as Epetra_Vector::Map() only returns an Epetra_BlockMap which cannot be
  // cast into an Epetra_Map, this workaround is needed.
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  return mesh_->getComplexNonOverlapMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
ModelEvaluator::
get_f_map() const
{
#ifdef _DEBUG_
    TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
    return mesh_->getComplexNonOverlapMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::
get_x_init () const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !x_.is_null() );
#endif
  return x_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::
get_p_init ( int l ) const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( 0, l );
#endif
  return p_init_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
ModelEvaluator::
get_p_map( int l ) const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( 0, l );
#endif
  return p_map_;
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
ModelEvaluator::
get_p_names( int l ) const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( 0, l );
#endif
  return p_names_;
}
// =============================================================================
Teuchos::RCP<Epetra_Operator>
ModelEvaluator::
create_W() const
{
  return Teuchos::rcp( new Ginla::JacobianOperator( mesh_, thickness_, keoFactory_ ) );
}
// =============================================================================
Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
ModelEvaluator::
create_WPrec() const
{
  Teuchos::RCP<Epetra_Operator> keoPrec =
          Teuchos::rcp( new Ginla::KeoRegularized( keoFactory_ ) );
  // bool is answer to: "Prec is already inverted?"
  // This needs to be set to TRUE to make sure that the constructor of
  //    NOX::Epetra::LinearSystemStratimikos
  // chooses a user-defined preconditioner.
  // Effectively, this boolean serves pretty well as a quirky switch for the
  // preconditioner if Piro is used.
  return Teuchos::rcp( new EpetraExt::ModelEvaluator::Preconditioner( keoPrec, true ) );
}
// ============================================================================
EpetraExt::ModelEvaluator::InArgs
ModelEvaluator::
createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  inArgs.setModelEvalDescription( "Ginzburg-Landau" );

  // We have *one* parameter vector with numParams_ parameters in it.
  inArgs.set_Np( 1 );

  inArgs.setSupports( IN_ARG_x, true );

  // for shifted matrix
  // TODO add support for operator shift
  inArgs.setSupports( IN_ARG_alpha, true );
  inArgs.setSupports( IN_ARG_beta, true );

  return inArgs;
}
// ============================================================================
EpetraExt::ModelEvaluator::OutArgs
ModelEvaluator::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  outArgs.setModelEvalDescription( "Ginzburg-Landau" );

  outArgs.set_Np_Ng( 1, 0 ); // one parameter vector, no objective function

  // support derivatives with respect to all parameters;
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
ModelEvaluator::
evalModel( const InArgs  & inArgs,
           const OutArgs & outArgs
         ) const
{
#ifdef GINLA_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm(*evalModelTime_);
#endif

  const double alpha = inArgs.get_alpha();
        double beta  = inArgs.get_beta();

  // From packages/piro/test/MockModelEval_A.cpp
  if (alpha==0.0 && beta==0.0)
  {
      //*out_ << "Ginla::ModelEvaluator Warning: alpha=beta=0 -- setting beta=1" << std::endl;
      beta = 1.0;
  }
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( alpha, 0.0 );
  TEUCHOS_ASSERT_EQUALITY( beta,  1.0 );
#endif

  const Teuchos::RCP<const Epetra_Vector> & x_in = inArgs.get_x();

  // get input arguments and make sure they are all right
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
#ifdef _DEBUG_
  TEUCHOS_ASSERT( !p_in.is_null() );
  for ( int k=0; k<p_in->MyLength(); k++ )
      TEUCHOS_ASSERT( !std::isnan( (*p_in)[k] ) );
#endif

  // p_current_ is used in getParameters.
  // Setting p_current_=p_in here is really a somewhat arbitrary choice.
  // The rationale is that the p-values which were last
  // used here are the "current" parameter values,
  // but there's actually no guarantee for it:
  // The evaluator could have been used for *anything.
  // Anyway, current_p_ is only used in this->getParameters().
  *p_current_ = *p_in;

  Teuchos::RCP<LOCA::ParameterVector> mvpParams = this->getParameters();

  const double temperature = (*p_in)[2];

  // compute F
  const Teuchos::RCP<Epetra_Vector> f_out = outArgs.get_f();
  if ( !f_out.is_null() )
  {
#ifdef GINLA_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor tm(*computeFTime_);
#endif
      this->computeF_( *x_in, mvpParams, temperature, *f_out );
  }

  // compute dF/dp
  const Teuchos::RCP<Epetra_MultiVector> dfdp_out =
      outArgs.get_DfDp(0).getMultiVector();
  if ( !dfdp_out.is_null() )
  {
#ifdef GINLA_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor tm(*computedFdpTime_);
#endif
     this->computedFdMu_( *x_in, mvpParams, *dfdp_out );
  }

  // fill jacobian
  const Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();
  if( !W_out.is_null() )
  {
#ifdef GINLA_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor tm(*fillJacobianTime_);
#endif
      Teuchos::RCP<Ginla::JacobianOperator> jac =
          Teuchos::rcp_dynamic_cast<Ginla::JacobianOperator>( W_out, true );
      jac->rebuild( mvpParams,
                    temperature,
                    x_in
                  );
  }

  // fill preconditioner
  const Teuchos::RCP<Epetra_Operator> WPrec_out = outArgs.get_WPrec();
  if( !WPrec_out.is_null() )
  {
#ifdef GINLA_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor tm(*fillPreconditionerTime_);
#endif
      Teuchos::RCP<Ginla::KeoRegularized> keoPrec =
          Teuchos::rcp_dynamic_cast<Ginla::KeoRegularized>( WPrec_out, true );
      keoPrec->rebuild( mvpParams );
  }

  return;
}
// ============================================================================
void
ModelEvaluator::
computeF_ ( const Epetra_Vector                             & x,
            const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams,
            const double                                      temperature,
            Epetra_Vector                                   & FVec
          ) const
{
  // build the KEO
  keoFactory_->updateParameters( mvpParams );
  const Teuchos::RCP<const Epetra_CrsMatrix> keoMatrix = keoFactory_->getKeo();

  // compute FVec = K*x
  TEUCHOS_ASSERT_EQUALITY( 0, keoMatrix->Apply( x, FVec ) );

  // add the nonlinear part (mass lumping)
#ifdef _DEBUG_
  TEUCHOS_ASSERT( FVec.Map().SameAs( x.Map() ) );
#endif

  const Epetra_Vector & controlVolumes = *(mesh_->getControlVolumes());

#ifdef _DEBUG_
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY( 2*controlVolumes.MyLength(), x.MyLength() );
#endif

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
      //     as suggested by mass lumping. This works if thickness(x_k) is available.
      //
      double alpha = controlVolumes[k] * (*thickness_)[k]
                     * ( (1.0-temperature) - x[2*k]*x[2*k] - x[2*k+1]*x[2*k+1] );
      // real part
      FVec[2*k]   += alpha * x[2*k];
      // imaginary part
      FVec[2*k+1] += alpha * x[2*k+1];
  }

  return;
}
// ============================================================================
void
ModelEvaluator::
computedFdMu_ ( const Epetra_Vector                             & x,
                const Teuchos::RCP<const LOCA::ParameterVector> & mvpParams,
                Epetra_MultiVector                              & FVec
              ) const
{
#ifdef _DEBUG_
  TEUCHOS_ASSERT_EQUALITY( 1, FVec.NumVectors() );
#endif
  // build the KEO
  keoFactory_->updateParameters( mvpParams );
  const Teuchos::RCP<const Epetra_CrsMatrix> dKdMuMatrix =
      keoFactory_->getKeoDMu();

  // compute FVec = K*x
  TEUCHOS_ASSERT_EQUALITY( 0, dKdMuMatrix->Apply( x, FVec ) );

  return;
}
// ============================================================================
Teuchos::RCP<Ginla::State>
ModelEvaluator::
createSavable( const Epetra_Vector & x ) const
{
    return Teuchos::rcp( new Ginla::State( x, mesh_ ) );
}
// =============================================================================
Teuchos::RCP<LOCA::ParameterVector>
ModelEvaluator::
getParameters() const
{
  // construct a LOCA::ParameterVector of the parameters
  Teuchos::RCP<LOCA::ParameterVector> p =
      Teuchos::rcp( new LOCA::ParameterVector() );

#ifdef _DEBUG_
  TEUCHOS_ASSERT( !p_names_.is_null() );
  TEUCHOS_ASSERT( !p_current_.is_null() );
#endif
  for ( int k=0; k<numParams_; k++ )
      p->addParameter( (*p_names_)[k], (*p_current_)[k] );

  return p;
}
// =============================================================================
} // namespace Ginla
