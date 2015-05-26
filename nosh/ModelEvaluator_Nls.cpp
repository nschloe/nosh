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

#include "nosh/ModelEvaluator_Nls.hpp"

#include "nosh/ScalarField_Virtual.hpp"
#include "nosh/ParameterMatrix_Virtual.hpp"
#include "nosh/JacobianOperator.hpp"
#include "nosh/KeoRegularized.hpp"
#include "nosh/StkMesh.hpp"

#include <string>
#include <map>

#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif

#include <Teuchos_VerboseObject.hpp>

namespace Nosh
{
namespace ModelEvaluator
{
// ============================================================================
Nls::
Nls(
  const Teuchos::RCP<const Nosh::StkMesh> &mesh,
  const Teuchos::RCP<Nosh::ParameterMatrix::Virtual> &keo,
  const Teuchos::RCP<Nosh::ParameterMatrix::Virtual> &dKeoDP,
  const Teuchos::RCP<const Nosh::ScalarField::Virtual> &scalarPotential,
  const double g,
  const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness,
  const Teuchos::RCP<const Epetra_Vector> &initialX
) :
  mesh_(mesh),
  scalarPotential_(scalarPotential),
  thickness_(thickness),
  x_init_(initialX),
  keo_(keo),
  dKeoDP_(dKeoDP),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  evalModelTime_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: Nls::evalModel"
        )),
  computeFTime_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: Nls::evalModel:compute F"
        )),
  computedFdpTime_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: Nls::evalModel:compute dF/dp"
        )),
  fillJacobianTime_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: Nls::evalModel:fill Jacobian"
        )),
  fillPreconditionerTime_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: Nls::evalModel::fill preconditioner"
        )),
#endif
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  p_map_(Teuchos::null),
  p_init_(Teuchos::null),
  p_names_(Teuchos::null)
{
  // Merge all of the parameters together.
  std::map<std::string, double> params;
  params["g"] = g;

  // This merges and discards new values if their keys are already in the list.
  std::map<std::string, double> spParams =
    scalarPotential_->getParameters();
  params.insert(spParams.begin(), spParams.end());

  // This merges and discards new values if their keys are already in the list.
  std::map<std::string, double> mbParams =
    keo_->getParameters();
  params.insert(mbParams.begin(), mbParams.end());

  // Out of this now complete list, create the entities that the EpetraExt::
  // Modelevaluator needs.
  const int numParams = params.size();
  p_map_ = Teuchos::rcp(new Epetra_LocalMap(numParams, 0, x_init_->Comm()));
  p_init_ = Teuchos::rcp(new Epetra_Vector(*p_map_));
  p_names_ = Teuchos::rcp(new Teuchos::Array<std::string>(numParams));
  int k = 0;
  for (auto it = params.begin(); it != params.end(); ++it) {
    (*p_names_)[k] = it->first;
    (*p_init_)[k] = it->second;
    k++;
  }

  return;
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
  TEUCHOS_ASSERT(!mesh_.is_null());
#endif
  return mesh_->getComplexNonOverlapMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Nls::
get_f_map() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(!mesh_.is_null());
#endif
  return mesh_->getComplexNonOverlapMap();
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Nls::
get_x_init() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(!x_init_.is_null());
#endif
  return x_init_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Nls::
get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      l != 0,
      "LOCA can only deal with one parameter vector."
      );
  return p_init_;
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Nls::
get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      l != 0,
      "LOCA can only deal with one parameter vector."
      );
  return p_map_;
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
Nls::
get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      l != 0,
      "LOCA can only deal with one parameter vector."
      );
  return p_names_;
}
// =============================================================================
Teuchos::RCP<Epetra_Operator>
Nls::
create_W() const
{
  return Teuchos::rcp(
      new Nosh::JacobianOperator(
        mesh_,
        scalarPotential_,
        thickness_,
        keo_
        )
      );
}
// =============================================================================
Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Nls::
create_WPrec() const
{
  Teuchos::RCP<Epetra_Operator> keoPrec = Teuchos::rcp(
      new Nosh::KeoRegularized(
        mesh_,
        thickness_,
        keo_
        )
      );
  // bool is answer to: "Prec is already inverted?"
  // This needs to be set to TRUE to make sure that the constructor of
  //    NOX::Epetra::LinearSystemStratimikos
  // chooses a user-defined preconditioner.
  // Effectively, this boolean serves pretty well as a quirky switch for the
  // preconditioner if Piro is used.
  return Teuchos::rcp(
      new EpetraExt::ModelEvaluator::Preconditioner(keoPrec, true)
      );
}
// ============================================================================
EpetraExt::ModelEvaluator::InArgs
Nls::
createInArgs() const
{
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;

  inArgs.setModelEvalDescription("Nonlinear Schr\"odinger");

  // We have *one* parameter vector with numParams_ parameters in it.
  inArgs.set_Np(1);

  inArgs.setSupports(IN_ARG_x, true);

  // for shifted matrix
  // TODO add support for operator shift
  inArgs.setSupports(IN_ARG_alpha, true);
  inArgs.setSupports(IN_ARG_beta, true);

  return inArgs;
}
// ============================================================================
EpetraExt::ModelEvaluator::OutArgs
Nls::
createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;

  outArgs.setModelEvalDescription("Nonlinear Schr\"odinger");

  outArgs.set_Np_Ng(1, 0); // one parameter vector, no objective function

  // support derivatives with respect to all parameters;
  outArgs.setSupports(OUT_ARG_DfDp,
                      0,
                      DerivativeSupport(DERIV_MV_BY_COL)
                     );

  outArgs.setSupports(OUT_ARG_f, true);
  outArgs.setSupports(OUT_ARG_W, true);
  outArgs.set_W_properties(
      DerivativeProperties(
        DERIV_LINEARITY_UNKNOWN, // DERIV_LINEARITY_NONCONST
        DERIV_RANK_DEFICIENT, // DERIV_RANK_FULL, DERIV_RANK_DEFICIENT
        false // supportsAdjoint
        )
      );

  outArgs.setSupports(OUT_ARG_WPrec, true);
  outArgs.set_WPrec_properties(
      DerivativeProperties(
        DERIV_LINEARITY_UNKNOWN,
        DERIV_RANK_FULL,
        false
        )
      );

  return outArgs;
}
// ============================================================================
void
Nls::
evalModel(
    const InArgs &inArgs,
    const OutArgs &outArgs
    ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm0(*evalModelTime_);
#endif

  const double alpha = inArgs.get_alpha();
  double beta = inArgs.get_beta();

  // From packages/piro/test/MockModelEval_A.cpp
  if (alpha == 0.0 && beta == 0.0) {
    beta = 1.0;
  }
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha, 0.0);
  TEUCHOS_ASSERT_EQUALITY(beta,  1.0);
#endif

  const Teuchos::RCP<const Epetra_Vector> &x_in = inArgs.get_x();
#ifndef NDEBUG
  TEUCHOS_ASSERT(!x_in.is_null());
#endif

  // Dissect inArgs.get_p(0) into parameter sublists.
  // Keep this in sync with get_p_init() where the splitting
  // is defined.
  const Teuchos::RCP<const Epetra_Vector> &p_in = inArgs.get_p(0);
#ifndef NDEBUG
  // Make sure the paremters aren't NaNs.
  TEUCHOS_ASSERT(!p_in.is_null());
  for (int k = 0; k < p_in->MyLength(); k++) {
    TEUCHOS_ASSERT(!std::isnan((*p_in)[k]));
  }
#endif

  // Fill the parameters into a std::map.
  const Teuchos::RCP<const Teuchos::Array<std::string> > paramNames =
    this->get_p_names(0);
  std::map<std::string, double> params;
  for (int k = 0; k < p_in->MyLength(); k++) {
    params[(*paramNames)[k]] = (*p_in)[k];
    //std::cout << (*paramNames)[k] << " " << (*p_in)[k] << std::endl;
  }


  // compute F
  const Teuchos::RCP<Epetra_Vector> &f_out = outArgs.get_f();
  if (!f_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm1(*computeFTime_);
#endif
    this->computeF_(*x_in, params, *f_out);
  }

  // Compute df/dp.
  const EpetraExt::ModelEvaluator::DerivativeMultiVector &derivMv =
    outArgs.get_DfDp(0).getDerivativeMultiVector();
  const Teuchos::RCP<Epetra_MultiVector> &dfdp_out =
    derivMv.getMultiVector();
  if (!dfdp_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm2(*computedFdpTime_);
#endif
    const Teuchos::Array<int> &paramIndices = derivMv.getParamIndexes();
    const int numDerivs = paramIndices.size();
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(numDerivs, dfdp_out->NumVectors());
#endif
    // Compute all derivatives.
    for (int k = 0; k < numDerivs; k++) {
      this->computeDFDP_(
          *x_in,
          params,
          (*paramNames)[paramIndices[k]],
          *(*dfdp_out)(k)
          );
    }
  }

  // Fill Jacobian.
  const Teuchos::RCP<Epetra_Operator> & W_out = outArgs.get_W();
  if(!W_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm3(*fillJacobianTime_);
#endif
    const Teuchos::RCP<Nosh::JacobianOperator> & jac =
      Teuchos::rcp_dynamic_cast<Nosh::JacobianOperator>(W_out, true);
    jac->rebuild(params, x_in);
  }

  // Fill preconditioner.
  const Teuchos::RCP<Epetra_Operator> & WPrec_out = outArgs.get_WPrec();
  if(!WPrec_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm4(*fillPreconditionerTime_);
#endif
    const Teuchos::RCP<Nosh::KeoRegularized> & keoPrec =
      Teuchos::rcp_dynamic_cast<Nosh::KeoRegularized>(WPrec_out, true);
    keoPrec->rebuild(params, *x_in);
  }

  return;
}
// ============================================================================
void
Nls::
computeF_(
    const Epetra_Vector &x,
    const std::map<std::string, double> & params,
    Epetra_Vector &FVec
    ) const
{
  // Compute FVec = K*x.
  keo_->setParameters(params);
  keo_->Apply(x, FVec);

  // Add the nonlinear part (mass lumping).
#ifndef NDEBUG
  TEUCHOS_ASSERT(FVec.Map().SameAs(x.Map()));
  TEUCHOS_ASSERT(!mesh_.is_null());
  TEUCHOS_ASSERT(!scalarPotential_.is_null());
  TEUCHOS_ASSERT(!thickness_.is_null());
#endif

  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());
  const int numMyPoints = controlVolumes.MyLength();

#ifndef NDEBUG
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, x.MyLength());
#endif

  const double g = params.at("g");

  const Epetra_Vector thicknessValues = thickness_->getV(params);
#ifndef NDEBUG
  TEUCHOS_ASSERT(controlVolumes.Map().SameAs(thicknessValues.Map()));
#endif

  const Epetra_Vector scalarPotentialValues = scalarPotential_->getV(params);
#ifndef NDEBUG
  TEUCHOS_ASSERT(controlVolumes.Map().SameAs(scalarPotentialValues.Map()));
#endif

  for (int k = 0; k < numMyPoints; k++) {
    // In principle, mass lumping here suggests to take
    //
    //   \int_{control volume}  thickness * f(psi).
    //
    // with f(psi) = psi[k] * ((1.0-T) - std::norm(psi[k])).
    // (a) A possible approximation for this is
    //
    //        |control volume| * average(thicknesses) * f(psi(x_k)).
    //
    //     This is loosely derived from the midpoint quadrature rule for
    //     triangles, i.e.,
    //
    //      \int_{triangle} f(x) ~=
    //        |triangle| * \sum_{edge midpoint} 1/3 * f(midpoint).
    //
    //     so all the values at the midpoints have the same weight (independent
    //     of whether the edge is long or short). This is then "generalized to
    //
    //      \int_{triangle} f(x)*a(x) ~=
    //        |triangle| *  \sum_{edge midpoint} 1/3 * f(midpoint)*a(midpoint),
    //
    //     or, as f(midpoint) is not available,
    //
    //      \int_{triangle} f(x)*a(x) ~=
    //        |triangle| * f(center of gravity)
    //          \sum_{edge midpoint} 1/3 * a(midpoint).
    //
    //     For general polynomals, this is then the above expression.
    //     Hence, do the equivalent of
    //
    //       res[k] += controlVolumes[k] * average(thicknesses)
    //               * psi[k] * (V + std::norm(psi[k]));
    //
    // (b) Another possible approximation is
    //
    //        |control volume| * thickness(x_k) * f(psi(x_k))
    //
    //     as suggested by mass lumping. This works if thickness(x_k) is
    //     available.
    //
    // The indexing here assumes that the local index K of controlVolume's map
    // is known to be local by thickness and scalarPotential and known to be
    // associated with that map.
    const double alpha =
      controlVolumes[k] * thicknessValues[k]* (
          scalarPotentialValues[k] + g * (x[2*k]*x[2*k] + x[2*k+1]*x[2*k+1])
          );
    // real and imaginary part
    FVec[2*k]   += alpha * x[2*k];
    FVec[2*k+1] += alpha * x[2*k+1];
  }

  return;
}
// ============================================================================
void
Nls::
computeDFDP_(
    const Epetra_Vector &x,
    const std::map<std::string, double> & params,
    const std::string & paramName,
    Epetra_Vector &FVec
    ) const
{
  // FVec = dK/dp * x.
  dKeoDP_->setParameters(params);
  dKeoDP_->Apply(x, FVec);

#ifndef NDEBUG
  TEUCHOS_ASSERT(FVec.Map().SameAs(x.Map()));
  TEUCHOS_ASSERT(!mesh_.is_null());
  TEUCHOS_ASSERT(!thickness_.is_null());
#endif
  const Epetra_Vector &controlVolumes = *(mesh_->getControlVolumes());

#ifndef NDEBUG
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY(2*controlVolumes.MyLength(), x.MyLength());
#endif

  const Epetra_Vector thicknessValues = thickness_->getV(params);
#ifndef NDEBUG
  TEUCHOS_ASSERT(controlVolumes.Map().SameAs(thicknessValues.Map()));
#endif

  if (paramName.compare("g") == 0) {
    for (int k = 0; k < controlVolumes.MyLength(); k++) {
      // This assumes that "g" is not a parameter in either of the
      // potentials.
      double alpha = controlVolumes[k] * thicknessValues[k]
                     * (x[2*k]*x[2*k] + x[2*k+1]*x[2*k+1]);
      // real and imaginary part
      FVec[2*k]   += alpha * x[2*k];
      FVec[2*k+1] += alpha * x[2*k+1];
    }
  } else {
    const Epetra_Vector scalarPotentialValues =
      scalarPotential_->getdVdP(params, paramName);
#ifndef NDEBUG
    TEUCHOS_ASSERT(controlVolumes.Map().SameAs(scalarPotentialValues.Map()));
#endif

    for (int k = 0; k < controlVolumes.MyLength(); k++) {
      const double alpha = controlVolumes[k] * thicknessValues[k]
                           * scalarPotentialValues[k];
      // real and imaginary part
      FVec[2*k]   += alpha * x[2*k];
      FVec[2*k+1] += alpha * x[2*k+1];
    }
  }

  return;
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
  for (int k = 0; k < numMyPoints; k++)
    res += controlVolumes[k] * (phi[2*k]*psi[2*k] + phi[2*k+1]*psi[2*k+1]);

  // Sum over all processors.
  double globalRes;
  TEUCHOS_ASSERT_EQUALITY(0, psi.Comm().SumAll(&res, &globalRes, 1));

  // normalize and return
  return globalRes / mesh_->getDomainVolume();
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
  for (int k = 0; k < numMyPoints; k++) {
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
}  // namespace ModelEvaluator
}  // namespace Nosh
