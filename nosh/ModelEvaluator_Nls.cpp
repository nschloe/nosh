// @HEADER
//
//    Nosh model evaluator.
//    Copyright (C) 2015  Nico Schlömer
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

#include "ModelEvaluator_Nls.hpp"

#include "ScalarField_Virtual.hpp"
#include "ParameterMatrix_Virtual.hpp"
#include "ParameterMatrix_Keo.hpp"
#include "ParameterMatrix_DKeoDP.hpp"
#include "JacobianOperator.hpp"
#include "KeoRegularized.hpp"
#include "StkMesh.hpp"
#include "Nosh_RealScalarProd.hpp"

#include <string>
#include <map>

#include <Thyra_DefaultPreconditioner.hpp>
#include <Thyra_ModelEvaluatorBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraThyraWrappers_decl.hpp>
#include <Thyra_TpetraLinearOp.hpp>
//#include <Thyra_MultiVectorAdapterBase.hpp>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_TimeMonitor.hpp>
#endif

#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

namespace Nosh
{
namespace ModelEvaluator
{
// ============================================================================
Nls::
Nls(
    const std::shared_ptr<const Nosh::StkMesh> &mesh,
    const std::shared_ptr<Nosh::VectorField::Virtual> &mvp,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &scalarPotential,
    const double g,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
    const std::shared_ptr<const Tpetra::Vector<double,int,int>> &initialX,
    const std::string & derivParameter
   ) :
  mesh_(mesh),
  mvp_(mvp),
  scalarPotential_(scalarPotential),
  thickness_(thickness),
  keo_(
      std::make_shared<Nosh::ParameterMatrix::Keo>(mesh_, thickness_, mvp_)
      ),
  dKeoDP_(
      std::make_shared<Nosh::ParameterMatrix::DKeoDP>(
        mesh_, thickness_, mvp_, derivParameter
        )
      ),
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
  p_names_(Teuchos::null),
  nominalValues_(this->createInArgs()),
  space_(createAlteredSpace())
{
  // Merge all of the parameters together.
  std::map<std::string, double> params;
  params["g"] = g;

  // This merges and discards new values if their keys are already in the list.
  auto spParams = scalarPotential_->getParameters();
  params.insert(spParams.begin(), spParams.end());

  // This merges and discards new values if their keys are already in the list.
  auto mbParams = keo_->getParameters();
  params.insert(mbParams.begin(), mbParams.end());

  // Out of this now complete list, create the entities that the Modelevaluator
  // needs.
  const int numParams = params.size();
  p_map_ = Teuchos::rcp(
      new Tpetra::Map<int,int>(
        numParams,
        0,
        Teuchos::rcp(mesh_->getComm())
        )
      );
  auto p_init = Thyra::createMember(this->get_p_space(0));
  p_names_ = Teuchos::rcp(new Teuchos::Array<std::string>(numParams));
  int k = 0;
  for (auto it = params.begin(); it != params.end(); ++it) {
    (*p_names_)[k] = it->first;
    Thyra::set_ele(k, it->second, p_init());
    k++;
  }

  // set nominal values
  auto xxx = Thyra::createConstVector(Teuchos::rcp(initialX), space_);
  nominalValues_.set_p(0, p_init);
  nominalValues_.set_x(xxx);

  return;
}
// ============================================================================
Nls::
~Nls()
{
}
// ============================================================================
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
Nls::
createAlteredSpace() const
{
  auto a = Thyra::createVectorSpace<double>(
      Teuchos::rcp(mesh_->getComplexNonOverlapMap()), true
      );
  // Use the Nosh scalar product. We still need to cast, cf.
  // <https://software.sandia.gov/bugzilla/show_bug.cgi?id=6355>.
  auto s = Teuchos::rcp_dynamic_cast<Thyra::ScalarProdVectorSpaceBase<double>>(a, true);
  auto sp = Teuchos::rcp(new Nosh::RealScalarProd<double>());
  //auto sp = Teuchos::rcp(new Thyra::EuclideanScalarProd<double>());
  s->setScalarProd(sp);

  return s;
}
// ============================================================================
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
Nls::
get_x_space() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(!space_.is_null());
#endif
  return space_;
}
// ============================================================================
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
Nls::
get_f_space() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(!space_.is_null());
#endif
  return space_;
}
// ============================================================================
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
Nls::
get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      l != 0,
      "LOCA can only deal with one parameter vector."
      );
  return Thyra::createVectorSpace<double>(p_map_);
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string>>
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
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
Nls::
get_g_space(int l) const
{
  (void) l;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      true,
      "Not implemented."
      );
  return Teuchos::null;
}
// =============================================================================
Thyra::ModelEvaluatorBase::InArgs<double>
Nls::
getNominalValues() const
{
  return nominalValues_;
}
// =============================================================================
Thyra::ModelEvaluatorBase::InArgs<double>
Nls::
getLowerBounds() const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      true,
      "Not implemented."
      );
  return this->createInArgs();
}
// =============================================================================
Thyra::ModelEvaluatorBase::InArgs<double>
Nls::
getUpperBounds() const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      true,
      "Not implemented."
      );
  return this->createInArgs();
}
// =============================================================================
Teuchos::RCP<Thyra::LinearOpBase<double>>
Nls::
create_W_op() const
{
  Teuchos::RCP<Tpetra::Operator<double,int,int>> jac = Teuchos::rcp(
        new Nosh::JacobianOperator(
          mesh_,
          scalarPotential_,
          thickness_,
          keo_
          )
        );

  return Thyra::createLinearOp(jac, space_, space_);
}
// =============================================================================
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
Nls::
get_W_factory() const
{
  Stratimikos::DefaultLinearSolverBuilder builder;

  auto p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "Belos");
  auto & belosList =
    p->sublist("Linear Solver Types")
    .sublist("Belos");
  //belosList.set("Solver Type", "MINRES");
  //belosList.set("Solver Type", "Pseudo Block GMRES");
  belosList.set("Solver Type", "Pseudo Block CG");

  auto & solverList =
    belosList.sublist("Solver Types")
    .sublist("Pseudo Block CG");
  solverList.set("Output Frequency", 1);
  solverList.set("Output Style", 1);
  solverList.set("Verbosity", 33);

  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  auto lowsFactory = builder.createLinearSolveStrategy("");

  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  return lowsFactory;
}
// =============================================================================
Teuchos::RCP<Thyra::PreconditionerBase<double>>
Nls::
create_W_prec() const
{
  const Teuchos::RCP<Tpetra::Operator<double,int,int>> keoPrec = Teuchos::rcp(
      new Nosh::KeoRegularized(
        mesh_,
        thickness_,
        mvp_
        )
      );
  auto keoT = Thyra::createLinearOp(keoPrec, space_, space_);
  return Thyra::nonconstUnspecifiedPrec(keoT);
}
// ============================================================================
void
Nls::
reportFinalPoint(
    const Thyra::ModelEvaluatorBase::InArgs<double> &finalPoint,
    const bool wasSolved
    )
{
  (void) finalPoint;
  (void) wasSolved;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      true,
      "Not implemented."
      );
}
// ============================================================================
Thyra::ModelEvaluatorBase::InArgs<double>
Nls::
createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<double> inArgs;

  inArgs.setModelEvalDescription("Nonlinear Schrödinger");

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
Thyra::ModelEvaluatorBase::OutArgs<double>
Nls::
createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<double> outArgs;

  outArgs.setModelEvalDescription("Nonlinear Schrödinger");

  outArgs.set_Np_Ng(1, 0); // one parameter vector, no objective function

  outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);

  // support derivatives with respect to all parameters;
  outArgs.setSupports(
      Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,
      0,
      DerivativeSupport(DERIV_MV_BY_COL)
      );

  outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
  //outArgs.set_W_properties(
  //    DerivativeProperties(
  //      DERIV_LINEARITY_UNKNOWN, // DERIV_LINEARITY_NONCONST
  //      DERIV_RANK_DEFICIENT, // DERIV_RANK_FULL, DERIV_RANK_DEFICIENT
  //      false // supportsAdjoint
  //      )
  //    );

  outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec);
  //outArgs.set_W_prec_properties(
  //    DerivativeProperties(
  //      DERIV_LINEARITY_UNKNOWN,
  //      DERIV_RANK_FULL,
  //      false
  //      )
  //    );

  return outArgs;
}
// ============================================================================
void
Nls::
evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<double> & inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<double> & outArgs
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

  const auto & x_in = inArgs.get_x();
#ifndef NDEBUG
  TEUCHOS_ASSERT(!x_in.is_null());
#endif
  // create corresponding tpetra vector
  auto x_in_tpetra =
    Thyra::TpetraOperatorVectorExtraction<double,int,int>::getConstTpetraVector(
        x_in
        );

  // Dissect inArgs.get_p(0) into parameter sublists.
  const auto & p_in = inArgs.get_p(0);
#ifndef NDEBUG
  TEUCHOS_ASSERT(!p_in.is_null());
#endif

#ifndef NDEBUG
  // Make sure the parameters aren't NaNs.
  for (int k = 0; k < p_in->space()->dim(); k++) {
    TEUCHOS_ASSERT(!std::isnan(Thyra::get_ele(*p_in, k)));
  }
#endif

  // Fill the parameters into a std::map.
  const auto paramNames = this->get_p_names(0);
  std::map<std::string, double> params;
  for (int k = 0; k < p_in->space()->dim(); k++) {
    params[(*paramNames)[k]] = Thyra::get_ele(*p_in, k);
    //std::cout << (*paramNames)[k] << " " << (*p_in)[k] << std::endl;
  }

  // compute F
  const auto & f_out = outArgs.get_f();
  if (!f_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm1(*computeFTime_);
#endif

    auto f_out_tpetra =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraVector(
          f_out
          );
    this->computeF_(
        *x_in_tpetra,
        params,
        *f_out_tpetra
        );
  }

  // Compute df/dp.
  const auto & derivMv = outArgs.get_DfDp(0).getDerivativeMultiVector();
  const auto & dfdp_out = derivMv.getMultiVector();
  if (!dfdp_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm2(*computedFdpTime_);
#endif
    auto dfdp_out_tpetra =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraMultiVector(
          dfdp_out
          );

    const int numAllParams = this->get_p_space(0)->dim();
    TEUCHOS_ASSERT_EQUALITY(
        numAllParams,
        dfdp_out_tpetra->getNumVectors()
        );
    // Compute all derivatives.
    for (int k = 0; k < numAllParams; k++) {
      this->computeDFDP_(
          *x_in_tpetra,
          params,
          (*paramNames)[k],
          *dfdp_out_tpetra->getVectorNonConst(k)
          );
    }
  }

  // Fill Jacobian.
  const auto & W_out = outArgs.get_W_op();
  if(!W_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm3(*fillJacobianTime_);
#endif
    auto W_outT =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraOperator(
          W_out
          );
    const auto & jac =
      Teuchos::rcp_dynamic_cast<Nosh::JacobianOperator>(W_outT, true);
    jac->rebuild(params, *x_in_tpetra);
  }

  // Fill preconditioner.
  const auto & WPrec_out = outArgs.get_W_prec();
  if(!WPrec_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm4(*fillPreconditionerTime_);
#endif
    auto WPrec_outT =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraOperator(
          WPrec_out->getNonconstUnspecifiedPrecOp()
          );
    const auto & keoPrec =
      Teuchos::rcp_dynamic_cast<Nosh::KeoRegularized>(WPrec_outT, true);
    keoPrec->rebuild(
        params,
        *x_in_tpetra
        );
  }

  return;
}
// ============================================================================
void
Nls::
computeF_(
    const Tpetra::Vector<double,int,int> &x,
    const std::map<std::string, double> & params,
    Tpetra::Vector<double,int,int> &FVec
    ) const
{
  // Compute FVec = K*x.
  keo_->setParameters(params);
  keo_->apply(x, FVec);

  auto xData = x.getData();
  auto fData = FVec.getDataNonConst();

  // Add the nonlinear part (mass lumping).
#ifndef NDEBUG
  TEUCHOS_ASSERT(FVec.getMap()->isSameAs(*x.getMap()));
  TEUCHOS_ASSERT(mesh_);
  TEUCHOS_ASSERT(scalarPotential_);
  TEUCHOS_ASSERT(thickness_);
#endif

  const auto & controlVolumes = *(mesh_->getControlVolumes());
  auto cData = controlVolumes.getData();

  const size_t numMyPoints = cData.size();

#ifndef NDEBUG
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, xData.size());
#endif

  const double g = params.at("g");

  const auto thicknessValues = thickness_->getV(params);
  auto tData = thicknessValues.getData();
#ifndef NDEBUG
  TEUCHOS_ASSERT(controlVolumes.getMap()->isSameAs(*thicknessValues.getMap()));
#endif

  const auto scalarPotentialValues = scalarPotential_->getV(params);
  auto sData = scalarPotentialValues.getData();
#ifndef NDEBUG
  TEUCHOS_ASSERT(controlVolumes.getMap()->isSameAs(*scalarPotentialValues.getMap()));
#endif

  for (size_t k = 0; k < numMyPoints; k++) {
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
      cData[k] * tData[k]* (
          sData[k] + g * (xData[2*k]*xData[2*k] + xData[2*k+1]*xData[2*k+1])
          );
    // real and imaginary part
    fData[2*k]   += alpha * xData[2*k];
    fData[2*k+1] += alpha * xData[2*k+1];
  }

  return;
}
// ============================================================================
void
Nls::
computeDFDP_(
    const Tpetra::Vector<double,int,int> &x,
    const std::map<std::string, double> & params,
    const std::string & paramName,
    Tpetra::Vector<double,int,int> &FVec
    ) const
{
  // FVec = dK/dp * x.
  dKeoDP_->setParameters(params);
  dKeoDP_->apply(x, FVec);

  auto xData = x.getData();
  auto fData = FVec.getDataNonConst();

#ifndef NDEBUG
  TEUCHOS_ASSERT(FVec.getMap()->isSameAs(*x.getMap()));
  TEUCHOS_ASSERT(mesh_);
  TEUCHOS_ASSERT(thickness_);
#endif
  const auto & controlVolumes = *(mesh_->getControlVolumes());
  auto cData = controlVolumes.getData();

#ifndef NDEBUG
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY(2*cData.size(), xData.size());
#endif

  const auto thicknessValues = thickness_->getV(params);
  auto tData = thicknessValues.getData();
#ifndef NDEBUG
  TEUCHOS_ASSERT(controlVolumes.getMap()->isSameAs(*thicknessValues.getMap()));
#endif

  if (paramName.compare("g") == 0) {
    for (int k = 0; k < cData.size(); k++) {
      // This assumes that "g" is not a parameter in either of the
      // potentials.
      double alpha = cData[k] * tData[k] *
        (xData[2*k]*xData[2*k] + xData[2*k+1]*xData[2*k+1]);
      // real and imaginary part
      fData[2*k]   += alpha * xData[2*k];
      fData[2*k+1] += alpha * xData[2*k+1];
    }
  } else {
    const auto scalarPotentialValues =
      scalarPotential_->getdVdP(params, paramName);
    auto sData = scalarPotentialValues.getData();
#ifndef NDEBUG
    TEUCHOS_ASSERT(controlVolumes.getMap()->isSameAs(
          *scalarPotentialValues.getMap()
          )
        );
#endif

    for (int k = 0; k < cData.size(); k++) {
      const double alpha = cData[k] * tData[k] * sData[k];
      // real and imaginary part
      fData[2*k]   += alpha * xData[2*k];
      fData[2*k+1] += alpha * xData[2*k+1];
    }
  }

  return;
}
// =============================================================================
double
Nls::
innerProduct(
    const Thyra::VectorBase<double> &phi,
    const Thyra::VectorBase<double> &psi
    ) const
{
  const auto & controlVolumes = *mesh_->getControlVolumes();

  auto cData = controlVolumes.getData();

  size_t numMyPoints = cData.size();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, phi.space()->dim());
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, psi.space()->dim());
#endif

  double res = 0.0;
  for (size_t k = 0; k < numMyPoints; k++) {
    double rePhi = Thyra::get_ele(phi, 2*k);
    double imPhi = Thyra::get_ele(phi, 2*k+1);
    double rePsi = Thyra::get_ele(psi, 2*k);
    double imPsi = Thyra::get_ele(psi, 2*k+1);
    res += cData[k] * (rePhi*rePsi + imPhi*imPsi);
  }

  // TODO
  //// Sum over all processors.
  //double globalRes;
  //TEUCHOS_ASSERT_EQUALITY(0,
  //    mesh_->getComm().reduceAll(
  //      &res,
  //      &globalRes,
  //      1
  //      )
  //    );
  //
  //// normalize and return
  //return globalRes / mesh_->getDomainVolume();

  return 0.0;
}
// =============================================================================
double
Nls::
gibbsEnergy(const Thyra::VectorBase<double> &psi) const
{
  const auto & controlVolumes = *mesh_->getControlVolumes();

  auto cData = controlVolumes.getData();

  size_t numMyPoints = cData.size();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, psi.space()->dim());
#endif

  double myEnergy = 0.0;
  for (size_t k = 0; k < numMyPoints; k++) {
    double re = Thyra::get_ele(psi, 2*k);
    double im = Thyra::get_ele(psi, 2*k+1);
    double alpha = re*re + im*im;
    myEnergy -= cData[k] * alpha * alpha;
  }

  // TODO
  //// Sum over all processors.
  //double globalEnergy;
  //TEUCHOS_ASSERT_EQUALITY(0, mesh_->getComm().SumAll(&myEnergy, &globalEnergy, 1));

  //// normalize and return
  //return globalEnergy / mesh_->getDomainVolume();

  return 0.0;
}
// =============================================================================
const std::shared_ptr<const Nosh::StkMesh>
Nls::
getMesh() const
{
  return mesh_;
}
// =============================================================================
}  // namespace ModelEvaluator
}  // namespace Nosh
