// @HEADER
//
//    Nosh model evaluator.
//    Copyright (C) 2010--2012  Nico Schl\"omer
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

#include "Nosh_ModelEvaluator_Nls.hpp"

#include "Nosh_ScalarField_Virtual.hpp"
#include "Nosh_MatrixBuilder_Virtual.hpp"
#include "Nosh_MatrixBuilder_Keo.hpp"
#include "Nosh_JacobianOperator.hpp"
#include "Nosh_KeoRegularized.hpp"
#include "Nosh_StkMesh.hpp"

#include <string>

#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  #include <Teuchos_TimeMonitor.hpp>
#endif

#include <Teuchos_VerboseObject.hpp>

namespace Nosh {
namespace ModelEvaluator {
// ============================================================================
Nls::
Nls(
  const Teuchos::RCP<const Nosh::StkMesh> &mesh,
  const Teuchos::RCP<const Nosh::MatrixBuilder::Virtual> &matrixBuilder,
  const Teuchos::RCP<const Nosh::ScalarField::Virtual> &scalarPotential,
  const double g,
  const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness,
  const Teuchos::RCP<const Epetra_Vector> &initialX
  ) :
  mesh_( mesh ),
  initial_g_( g ),
  scalarPotential_( scalarPotential ),
  thickness_( thickness ),
  x_init_( initialX ),
  p_latest_(Teuchos::null),
  matrixBuilder_(matrixBuilder),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  evalModelTime_( Teuchos::TimeMonitor::getNewTimer(
                    "Nosh: Nls::evalModel" ) ),
  computeFTime_( Teuchos::TimeMonitor::getNewTimer(
                   "Nosh: Nls::evalModel:compute F" ) ),
  computedFdpTime_( Teuchos::TimeMonitor::getNewTimer(
                      "Nosh: Nls::evalModel:compute dF/dp" ) ),
  fillJacobianTime_( Teuchos::TimeMonitor::getNewTimer(
                       "Nosh: Nls::evalModel:fill Jacobian" ) ),
  fillPreconditionerTime_( Teuchos::TimeMonitor::getNewTimer(
                             "Nosh: Nls::evalModel::fill preconditioner" ) ),
#endif
  out_( Teuchos::VerboseObjectBase::getDefaultOStream() )
{
}
// ============================================================================
Nls::
~Nls()
{
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Nls::
get_x_map() const
{
  // It is a bit of an assumption that x_ actually has this map, but
  // as Epetra_Vector::Map() only returns an Epetra_BlockMap which cannot be
  // cast into an Epetra_Map, this workaround is needed.
#ifndef NDEBUG
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  return mesh_->getComplexNonOverlapMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Nls::
get_f_map() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !mesh_.is_null() );
#endif
  return mesh_->getComplexNonOverlapMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Nls::
get_x_init() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( !x_init_.is_null() );
#endif
  return x_init_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Nls::
get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(l != 0,
                              "LOCA can only deal with one parameter vector.");
  // Put all of the parameters in one vector and remember where the
  // values are put. This information is used later on to distribute
  // incoming parameter vectors into where the parameters are
  // actually stored.
  Teuchos::RCP<Epetra_Vector> p_init =
    Teuchos::rcp(new Epetra_Vector(*this->get_p_map(l)));
  int k = 0;
  Teuchos::RCP<const Teuchos::Array<double> > p;
  // Local parameters:
  (*p_init)[k++] = initial_g_;
  // Scalar potential parameters:
  p = scalarPotential_->get_p_init();
  for (int i=0; i<p->length(); i++)
    (*p_init)[k++] = (*p)[i];
  // Energy operator parameters:
  p = matrixBuilder_->get_p_init();
  for (int i=0; i<p->length(); i++)
    (*p_init)[k++] = (*p)[i];

  return p_init;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Nls::
get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(l != 0,
                              "LOCA can only deal with one parameter vector.");
  int totalNumParams = 1 // local parameters
                     + scalarPotential_->get_p_init()->length() // scalar potential
                     + matrixBuilder_->get_p_init()->length(); // posdef operator

  return Teuchos::rcp(new Epetra_LocalMap(totalNumParams, 0, x_init_->Comm()));
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
Nls::
get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(l != 0,
                              "LOCA can only deal with one parameter vector.");
  int totalNumParams = 1 // local parameters
                     + scalarPotential_->get_p_names()->length() // scalar potential
                     + matrixBuilder_->get_p_names()->length(); // posdef operator

  Teuchos::RCP<Teuchos::Array<std::string> > p_names =
    Teuchos::rcp(new Teuchos::Array<std::string>(totalNumParams));
  int k = 0;
  Teuchos::RCP<const Teuchos::Array<std::string> > p;
  // Local parameters:
  (*p_names)[k++] = "g";
  // Scalar potential parameters:
  p = scalarPotential_->get_p_names();
  for (int i=0; i<p->length(); i++)
    (*p_names)[k++] = (*p)[i];
  // Vector potential parameters:
  p = matrixBuilder_->get_p_names();
  for (int i=0; i<p->length(); i++)
    (*p_names)[k++] = (*p)[i];

  return p_names;
}
// =============================================================================
Teuchos::RCP<Epetra_Operator>
Nls::
create_W() const
{
  return Teuchos::rcp(new Nosh::JacobianOperator(mesh_,
                                                 scalarPotential_,
                                                 thickness_,
                                                 matrixBuilder_));
}
// =============================================================================
Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Nls::
create_WPrec() const
{
  Teuchos::RCP<Epetra_Operator> keoPrec =
    Teuchos::rcp(new Nosh::KeoRegularized(mesh_,
                                          thickness_,
                                          matrixBuilder_));
  // bool is answer to: "Prec is already inverted?"
  // This needs to be set to TRUE to make sure that the constructor of
  //    NOX::Epetra::LinearSystemStratimikos
  // chooses a user-defined preconditioner.
  // Effectively, this boolean serves pretty well as a quirky switch for the
  // preconditioner if Piro is used.
  return Teuchos::rcp(new EpetraExt::ModelEvaluator::Preconditioner(keoPrec,
                                                                    true));
                                                                    //false));
}
// ============================================================================
EpetraExt::ModelEvaluator::InArgs
Nls::
createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  inArgs.setModelEvalDescription( "Nonlinear Schrödinger" );

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
Nls::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  outArgs.setModelEvalDescription( "Nonlinear Schrödinger" );

  outArgs.set_Np_Ng( 1, 0 ); // one parameter vector, no objective function

  // support derivatives with respect to all parameters;
  outArgs.setSupports(OUT_ARG_DfDp,
                      0,
                      DerivativeSupport( DERIV_MV_BY_COL )
                      );

  outArgs.setSupports( OUT_ARG_f, true );
  outArgs.setSupports( OUT_ARG_W, true );
  outArgs.set_W_properties(DerivativeProperties(DERIV_LINEARITY_UNKNOWN, // DERIV_LINEARITY_NONCONST
                                                DERIV_RANK_DEFICIENT, // DERIV_RANK_FULL, DERIV_RANK_DEFICIENT
                                                false // supportsAdjoint
                                                )
                           );

  outArgs.setSupports(OUT_ARG_WPrec, true);
  outArgs.set_WPrec_properties(DerivativeProperties(DERIV_LINEARITY_UNKNOWN,
                                                    DERIV_RANK_FULL,
                                                    false
                                                    )
                               );

  return outArgs;
}
// ============================================================================
void
Nls::
evalModel(const InArgs &inArgs,
          const OutArgs &outArgs
          ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm0( *evalModelTime_ );
#endif

  const double alpha = inArgs.get_alpha();
  double beta = inArgs.get_beta();

  // From packages/piro/test/MockModelEval_A.cpp
  if (alpha==0.0 && beta==0.0)
    beta = 1.0;
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha, 0.0);
  TEUCHOS_ASSERT_EQUALITY(beta,  1.0);
#endif

  const Teuchos::RCP<const Epetra_Vector> &x_in = inArgs.get_x();
#ifndef NDEBUG
  TEUCHOS_ASSERT( !x_in.is_null() );
#endif

  // Dissect inArgs.get_p(0) into parameter sublists.
  // Keep this in sync with get_p_init() where the splitting
  // is defined.
  const Teuchos::RCP<const Epetra_Vector> &p_in = inArgs.get_p(0);
#ifndef NDEBUG
    TEUCHOS_ASSERT( !p_in.is_null() );
    for (int k=0; k<p_in->MyLength(); k++)
      TEUCHOS_ASSERT( !std::isnan( (*p_in)[k] ) );
#endif
  int i = 0;
  // Gather g.
  const double g = (*p_in)[i++];
  // Gather scalar potential parameters.
  const int numSpParams = scalarPotential_->get_p_init()->length();
  Teuchos::Array<double> spParams(numSpParams);
  for (int k=0; k<numSpParams; k++)
    spParams[k] = (*p_in)[i++];
  // Gather vector potential parameters.
  const int numMvpParams = matrixBuilder_->get_p_init()->length();
  Teuchos::Array<double> eoParams(numMvpParams);
  for (int k=0; k<numMvpParams; k++)
    eoParams[k] = (*p_in)[i++];
  // Make sure we arrived at the end of the vector.
  TEUCHOS_ASSERT_EQUALITY(i, p_in->MyLength());

  // Store "current" parameters, used in this->getParameters().
  // Setting p_latest__=p_in here is really a somewhat arbitrary choice.
  // The rationale is that the p-values which were last
  // used here are the "current" parameter values,
  // but there's actually no guarantee for it:
  // The evaluator could have been used for *anything*.
  // Anyways, it's only used in this->get_p_latest (for NOX::Observer).
  p_latest_ = Teuchos::rcp(new Epetra_Vector(*p_in));

  // compute F
  const Teuchos::RCP<Epetra_Vector> &f_out = outArgs.get_f();
  if ( !f_out.is_null() )
  {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm1( *computeFTime_ );
#endif
    this->computeF_(*x_in, g, spParams, eoParams, *f_out);
  }

  // Compute df/dp.
  const EpetraExt::ModelEvaluator::DerivativeMultiVector &derivMv =
    outArgs.get_DfDp(0).getDerivativeMultiVector();
  const Teuchos::RCP<Epetra_MultiVector> &dfdp_out =
    derivMv.getMultiVector();
  if ( !dfdp_out.is_null() )
  {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm2( *computedFdpTime_ );
#endif
    const Teuchos::Array<int> &paramIndices = derivMv.getParamIndexes();
    const int numDerivs = paramIndices.size();
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(numDerivs, dfdp_out->NumVectors());
#endif
    // Come up with something better here, e.g.,
    // Initialize array of function pointers for derivative getters
    // so they can easily be employed in a loop?
    // Check out
    // http://www.parashift.com/c++-faq-lite/pointers-to-members.html#faq-33.7
    for (int k=0; k<numDerivs; k++)
    {
      if (paramIndices[k] < 1)
        this->computeDFDg_(*x_in,
                           *(*dfdp_out)(k));
      else if (paramIndices[k] < 1 + numSpParams)
        this->computeDFDPpotential_(*x_in,
                                    spParams,
                                    paramIndices[k] - 1,
                                    *(*dfdp_out)(k));
      else if (paramIndices[k] < 1 + numSpParams + numMvpParams)
        this->computeDFDPeo_(*x_in,
                             eoParams,
                             paramIndices[k] - 1 - numSpParams,
                             *(*dfdp_out)(k));
      else
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                                    "Illegal parameter index " << paramIndices[k] << ". Abort.");
    }
  }

  // Fill Jacobian.
  const Teuchos::RCP<Epetra_Operator> & W_out = outArgs.get_W();
  if( !W_out.is_null() )
  {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm3( *fillJacobianTime_ );
#endif
    const Teuchos::RCP<Nosh::JacobianOperator> & jac =
      Teuchos::rcp_dynamic_cast<Nosh::JacobianOperator>(W_out, true);
    jac->rebuild(g, spParams, eoParams, x_in);
  }

  // Fill preconditioner.
  const Teuchos::RCP<Epetra_Operator> & WPrec_out = outArgs.get_WPrec();
  if( !WPrec_out.is_null() )
  {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm4( *fillPreconditionerTime_ );
#endif
    const Teuchos::RCP<Nosh::KeoRegularized> & keoPrec =
      Teuchos::rcp_dynamic_cast<Nosh::KeoRegularized>(WPrec_out, true);
    keoPrec->rebuild(g, eoParams, *x_in);
  }

  return;
}
// ============================================================================
void
Nls::
computeF_(const Epetra_Vector &x,
          const double g,
          const Teuchos::Array<double> & spParams,
          const Teuchos::Array<double> & eoParams,
          Epetra_Vector &FVec
          ) const
{
  // Compute FVec = K*x.
  matrixBuilder_->apply(eoParams, x, FVec);

  // Add the nonlinear part (mass lumping).
#ifndef NDEBUG
  TEUCHOS_ASSERT( FVec.Map().SameAs( x.Map() ) );
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !scalarPotential_.is_null() );
  TEUCHOS_ASSERT( !thickness_.is_null() );
#endif

  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());
  const int numMyPoints = controlVolumes.MyLength();

#ifndef NDEBUG
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, x.MyLength());
#endif

  for (int k=0; k<numMyPoints; k++)
  {
    // In principle, mass lumping here suggests to take
    //
    //   \int_{control volume}  thickness * f(psi).
    //
    // with f(psi) = psi[k] * ( (1.0-T) - std::norm(psi[k]) ).
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
    //       res[k] += controlVolumes[k] * average(thicknesses) * psi[k] * ( V + std::norm(psi[k]) );
    //
    // (b) Another possible approximation is
    //
    //        |control volume| * thickness(x_k) * f(psi(x_k))
    //
    //     as suggested by mass lumping. This works if thickness(x_k) is available.
    //
    // The indexing here assumes that the local index K of controlVolume's map
    // is known to be local by thickness and scalarPotential and known to be
    // associated with that map.
    double alpha = controlVolumes[k] * thickness_->getV(k)
                 * (scalarPotential_->getV(k, spParams) + g * (x[2*k]*x[2*k] + x[2*k+1]*x[2*k+1]));
    // real and imaginary part
    FVec[2*k]   += alpha * x[2*k];
    FVec[2*k+1] += alpha * x[2*k+1];
  }

  return;
}
// ============================================================================
void
Nls::
computeDFDg_(const Epetra_Vector &x,
             Epetra_Vector &FVec
             ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( FVec.Map().SameAs( x.Map() ) );
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !thickness_.is_null() );
#endif
  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());

#ifndef NDEBUG
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY( 2*controlVolumes.MyLength(), x.MyLength() );
#endif

  for ( int k=0; k<controlVolumes.MyLength(); k++ )
  {
    double alpha = controlVolumes[k] * thickness_->getV(k)
                 * (x[2*k]*x[2*k] + x[2*k+1]*x[2*k+1]);
    // real and imaginary part
    FVec[2*k]   = alpha * x[2*k];
    FVec[2*k+1] = alpha * x[2*k+1];
  }

  return;
}
// ============================================================================
void
Nls::
computeDFDPpotential_(const Epetra_Vector &x,
                      const Teuchos::Array<double> & spParams,
                      int paramIndex,
                      Epetra_Vector &FVec
                      ) const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT( FVec.Map().SameAs( x.Map() ) );
  TEUCHOS_ASSERT( !mesh_.is_null() );
  TEUCHOS_ASSERT( !thickness_.is_null() );
#endif

  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());

#ifndef NDEBUG
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY(2*controlVolumes.MyLength(), x.MyLength());
#endif

  for ( int k=0; k<controlVolumes.MyLength(); k++ )
  {
    double alpha = controlVolumes[k] * thickness_->getV(k)
                 * scalarPotential_->getdVdP(k, paramIndex, spParams);
    // real and imaginary part
    FVec[2*k]   = alpha * x[2*k];
    FVec[2*k+1] = alpha * x[2*k+1];
  }

  return;
}
// ============================================================================
void
Nls::
computeDFDPeo_(const Epetra_Vector &x,
                const Teuchos::Array<double> &eoParams,
                int paramIndex,
                Epetra_Vector &FVec
                ) const
{
  // FVec = dK/dp * x.
  matrixBuilder_->applyDKDp(eoParams, paramIndex, x, FVec);

  return;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Nls::
get_p_latest() const
{
  // This is fetcher routine to make sure that the NOX::Observer
  // can print the parameter values.
  return p_latest_;
}
// =============================================================================
double
Nls::
innerProduct(const Epetra_Vector &phi,
             const Epetra_Vector &psi
             ) const
{
  const Epetra_Vector &controlVolumes = *mesh_->getControlVolumes();

  int numMyPoints = controlVolumes.Map().NumMyPoints();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, phi.MyLength());
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, psi.MyLength());
#endif

  double res = 0.0;
  for (int k=0; k<numMyPoints; k++)
    res += controlVolumes[k] * (phi[2*k]*psi[2*k] + phi[2*k+1]*psi[2*k+1]);

  // Sum over all processors.
  double globalRes;
  TEUCHOS_ASSERT_EQUALITY(0, psi.Comm().SumAll(&res, &globalRes, 1));

  // normalize and return
  return globalRes / mesh_->getDomainVolume();
}
// ============================================================================
double
Nls::
normalizedScaledL2Norm(const Epetra_Vector &psi) const
{
  return sqrt(this->innerProduct(psi, psi));
}
// =============================================================================
double
Nls::
gibbsEnergy(const Epetra_Vector &psi) const
{
  const Epetra_Vector &controlVolumes = *mesh_->getControlVolumes();

  int numMyPoints = controlVolumes.Map().NumMyPoints();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, psi.MyLength());
#endif

  double myEnergy = 0.0;
  double alpha;
  for (int k=0; k<numMyPoints; k++)
  {
    alpha = psi[2*k]*psi[2*k] + psi[2*k+1]*psi[2*k+1];
    myEnergy -= controlVolumes[k] * alpha * alpha;
  }

  // Sum over all processors.
  double globalEnergy;
  TEUCHOS_ASSERT_EQUALITY(0, psi.Comm().SumAll(&myEnergy, &globalEnergy, 1));

  // normalize and return
  return globalEnergy / mesh_->getDomainVolume();
}
// =============================================================================
const Teuchos::RCP<const Nosh::StkMesh>
Nls::
getMesh() const
{
  return mesh_;
}
// =============================================================================
} // namespace ModelEvaluator
} // namespace Nosh
