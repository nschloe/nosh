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

#include "model_evaluator_nls.hpp"

#include "scalar_field_base.hpp"
#include "parameter_matrix_base.hpp"
#include "parameter_matrix_keo.hpp"
#include "parameter_matrix_dkeo_dp.hpp"
#include "jacobian_operator.hpp"
#include "keo_regularized.hpp"
#include "mesh.hpp"
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

namespace nosh
{
namespace model_evaluator
{
// ============================================================================
nls::
nls(
    const std::shared_ptr<const nosh::mesh> &_mesh,
    const std::shared_ptr<nosh::vector_field::base> &mvp,
    const std::shared_ptr<const nosh::scalar_field::base> &scalar_potential,
    const double g,
    const std::shared_ptr<const nosh::scalar_field::base> &thickness,
    const std::shared_ptr<const Tpetra::Vector<double,int,int>> &initial_x,
    const std::string & deriv_parameter
   ) :
  mesh_(_mesh),
  mvp_(mvp),
  scalar_potential_(scalar_potential),
  thickness_(thickness),
  keo_(
      std::make_shared<nosh::parameter_matrix::keo>(mesh_, thickness_, mvp_)
      ),
  dkeo_dp_(
      std::make_shared<nosh::parameter_matrix::DkeoDP>(
        mesh_, thickness_, mvp_, deriv_parameter
        )
      ),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  eval_model_time_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: nls::eval_model"
        )),
  compute_f_time_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: nls::eval_model:compute F"
        )),
  compute_dfdp_time_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: nls::eval_model:compute dF/dp"
        )),
  fill_jacobian_time_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: nls::eval_model:fill Jacobian"
        )),
  fill_preconditioner_time_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: nls::eval_model::fill preconditioner"
        )),
#endif
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  p_map_(Teuchos::null),
  p_names_(Teuchos::null),
  nominal_values_(this->createInArgs()),
  space_(createAlteredSpace())
{
  // Merge all of the parameters together.
  std::map<std::string, double> params;
  params["g"] = g;

  // This merges and discards new values if their keys are already in the list.
  auto spParams = scalar_potential_->get_scalar_parameters();
  params.insert(spParams.begin(), spParams.end());

  // This merges and discards new values if their keys are already in the list.
  auto mbParams = keo_->get_scalar_parameters();
  params.insert(mbParams.begin(), mbParams.end());

  // Out of this now complete list, create the entities that the Modelevaluator
  // needs.
  const int numParams = params.size();
  p_map_ = Teuchos::rcp(
      new Tpetra::Map<int,int>(
        numParams,
        0,
        Teuchos::rcp(mesh_->comm)
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
  auto xxx = Thyra::createConstVector(Teuchos::rcp(initial_x), space_);
  nominal_values_.set_p(0, p_init);
  nominal_values_.set_x(xxx);

  return;
}
// ============================================================================
nls::
~nls()
{
}
// ============================================================================
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
nls::
createAlteredSpace() const
{
  auto a = Thyra::createVectorSpace<double>(
      Teuchos::rcp(mesh_->complex_map())
      );
#if 0
  // Use the Nosh scalar product. We still need to cast, cf.
  // <https://software.sandia.gov/bugzilla/show_bug.cgi?id=6355>.
  auto s = Teuchos::rcp_dynamic_cast<Thyra::ScalarProdVectorSpaceBase<double>>(a, true);
  auto sp = Teuchos::rcp(new nosh::RealScalarProd<double>());
  //auto sp = Teuchos::rcp(new Thyra::EuclideanScalarProd<double>());
  s->setScalarProd(sp);

  return s;
#endif
  return a;
}
// ============================================================================
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
nls::
get_x_space() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(!space_.is_null());
#endif
  return space_;
}
// ============================================================================
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
nls::
get_f_space() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(!space_.is_null());
#endif
  return space_;
}
// ============================================================================
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
nls::
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
nls::
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
nls::
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
Teuchos::ArrayView<const std::string>
nls::
get_g_names(int j) const
{
  (void) j;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      true,
      "Not implemented."
      );
  return Teuchos::null;
}
// =============================================================================
Thyra::ModelEvaluatorBase::InArgs<double>
nls::
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
nls::
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
nls::
create_W_op() const
{
  Teuchos::RCP<Tpetra::Operator<double,int,int>> jac = Teuchos::rcp(
        new nosh::jacobian_operator(
          mesh_,
          scalar_potential_,
          thickness_,
          keo_
          )
        );

  return Thyra::createLinearOp(jac, space_, space_);
}
// =============================================================================
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
nls::
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
nls::
create_W_prec() const
{
  const Teuchos::RCP<Tpetra::Operator<double,int,int>> keoPrec = Teuchos::rcp(
      new nosh::keo_regularized(
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
nls::
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
nls::
createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<double> in_args;

  in_args.setModelEvalDescription("Nonlinear Schrödinger");

  // We have *one* parameter vector with numParams_ parameters in it.
  in_args.set_Np(1);

  in_args.setSupports(IN_ARG_x, true);

  // for shifted matrix
  // TODO add support for operator shift
  in_args.setSupports(IN_ARG_alpha, true);
  in_args.setSupports(IN_ARG_beta, true);

  return in_args;
}
// ============================================================================
Thyra::ModelEvaluatorBase::OutArgs<double>
nls::
createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<double> out_args;

  out_args.setModelEvalDescription("Nonlinear Schrödinger");

  out_args.set_Np_Ng(1, 0); // one parameter vector, no objective function

  out_args.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);

  // support derivatives with respect to all parameters;
  out_args.setSupports(
      Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,
      0,
      DerivativeSupport(DERIV_MV_BY_COL)
      );

  out_args.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
  //out_args.set_W_properties(
  //    DerivativeProperties(
  //      DERIV_LINEARITY_UNKNOWN, // DERIV_LINEARITY_NONCONST
  //      DERIV_RANK_DEFICIENT, // DERIV_RANK_FULL, DERIV_RANK_DEFICIENT
  //      false // supportsAdjoint
  //      )
  //    );

  out_args.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec);
  //out_args.set_W_prec_properties(
  //    DerivativeProperties(
  //      DERIV_LINEARITY_UNKNOWN,
  //      DERIV_RANK_FULL,
  //      false
  //      )
  //    );

  return out_args;
}
// ============================================================================
void
nls::
evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<double> & in_args,
    const Thyra::ModelEvaluatorBase::OutArgs<double> & out_args
    ) const
{
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor tm0(*eval_model_time_);
#endif

  const double alpha = in_args.get_alpha();
  double beta = in_args.get_beta();

  // From packages/piro/test/MockModelEval_A.cpp
  if (alpha == 0.0 && beta == 0.0) {
    beta = 1.0;
  }
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(alpha, 0.0);
  TEUCHOS_ASSERT_EQUALITY(beta,  1.0);
#endif

  const auto & x_in = in_args.get_x();
#ifndef NDEBUG
  TEUCHOS_ASSERT(!x_in.is_null());
#endif
  // create corresponding tpetra vector
  auto x_in_tpetra =
    Thyra::TpetraOperatorVectorExtraction<double,int,int>::getConstTpetraVector(
        x_in
        );

  // Dissect in_args.get_p(0) into parameter sublists.
  const auto & p_in = in_args.get_p(0);
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
  const auto param_names = this->get_p_names(0);
  std::map<std::string, double> params;
  for (int k = 0; k < p_in->space()->dim(); k++) {
    params[(*param_names)[k]] = Thyra::get_ele(*p_in, k);
  }

  // compute F
  const auto & f_out = out_args.get_f();
  if (!f_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm1(*compute_f_time_);
#endif

    auto f_out_tpetra =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraVector(
          f_out
          );
    this->compute_f_(
        *x_in_tpetra,
        params,
        *f_out_tpetra
        );
  }

  // Compute df/dp.
  const auto & derivMv = out_args.get_DfDp(0).getDerivativeMultiVector();
  const auto & dfdp_out = derivMv.getMultiVector();
  if (!dfdp_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm2(*compute_dfdp_time_);
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
          (*param_names)[k],
          *dfdp_out_tpetra->getVectorNonConst(k)
          );
    }
  }

  // Fill Jacobian.
  const auto & W_out = out_args.get_W_op();
  if(!W_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm3(*fill_jacobian_time_);
#endif
    auto W_outT =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraOperator(
          W_out
          );
    const auto & jac =
      Teuchos::rcp_dynamic_cast<nosh::jacobian_operator>(W_outT, true);
    jac->rebuild(params, *x_in_tpetra);
  }

  // Fill preconditioner.
  const auto & WPrec_out = out_args.get_W_prec();
  if(!WPrec_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm4(*fill_preconditioner_time_);
#endif
    auto WPrec_outT =
      Thyra::TpetraOperatorVectorExtraction<double,int,int>::getTpetraOperator(
          WPrec_out->getNonconstUnspecifiedPrecOp()
          );
    const auto & keoPrec =
      Teuchos::rcp_dynamic_cast<nosh::keo_regularized>(WPrec_outT, true);
    keoPrec->rebuild(
        params,
        *x_in_tpetra
        );
  }

  return;
}
// ============================================================================
void
nls::
compute_f_(
    const Tpetra::Vector<double,int,int> &x,
    const std::map<std::string, double> & params,
    Tpetra::Vector<double,int,int> &f_vec
    ) const
{
  // Compute f_vec = K*x.
  keo_->set_parameters(params, {});
  keo_->apply(x, f_vec);

  auto x_data = x.getData();
  auto f_data = f_vec.getDataNonConst();

  // Add the nonlinear part (mass lumping).
#ifndef NDEBUG
  TEUCHOS_ASSERT(f_vec.getMap()->isSameAs(*x.getMap()));
  TEUCHOS_ASSERT(mesh_);
  TEUCHOS_ASSERT(scalar_potential_);
  TEUCHOS_ASSERT(thickness_);
#endif

  const auto & control_volumes = *(mesh_->control_volumes());
  auto c_data = control_volumes.getData();

  const size_t num_my_points = c_data.size();

#ifndef NDEBUG
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY(2*num_my_points, x_data.size());
#endif

  const double g = params.at("g");

  const auto thickness_values = thickness_->get_v(params);
  auto t_data = thickness_values.getData();
#ifndef NDEBUG
  TEUCHOS_ASSERT(control_volumes.getMap()->isSameAs(*thickness_values.getMap()));
#endif

  const auto scalar_potential_values = scalar_potential_->get_v(params);
  auto s_data = scalar_potential_values.getData();
#ifndef NDEBUG
  TEUCHOS_ASSERT(control_volumes.getMap()->isSameAs(*scalar_potential_values.getMap()));
#endif

  for (size_t k = 0; k < num_my_points; k++) {
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
    //       res[k] += control_volumes[k] * average(thicknesses)
    //               * psi[k] * (V + std::norm(psi[k]));
    //
    // (b) Another possible approximation is
    //
    //        |control volume| * thickness(x_k) * f(psi(x_k))
    //
    //     as suggested by mass lumping. This works if thickness(x_k) is
    //     available.
    //
    // The indexing here assumes that the local index K of control_volume's map
    // is known to be local by thickness and scalar_potential and known to be
    // associated with that map.
    const double alpha =
      c_data[k] * t_data[k]* (
          s_data[k] + g * (x_data[2*k]*x_data[2*k] + x_data[2*k+1]*x_data[2*k+1])
          );
    // real and imaginary part
    f_data[2*k]   += alpha * x_data[2*k];
    f_data[2*k+1] += alpha * x_data[2*k+1];
  }

  return;
}
// ============================================================================
void
nls::
computeDFDP_(
    const Tpetra::Vector<double,int,int> &x,
    const std::map<std::string, double> & params,
    const std::string & param_name,
    Tpetra::Vector<double,int,int> &f_vec
    ) const
{
  // f_vec = dK/dp * x.
  dkeo_dp_->set_parameters(params, {});
  dkeo_dp_->apply(x, f_vec);

  auto x_data = x.getData();
  auto f_data = f_vec.getDataNonConst();

#ifndef NDEBUG
  TEUCHOS_ASSERT(f_vec.getMap()->isSameAs(*x.getMap()));
  TEUCHOS_ASSERT(mesh_);
  TEUCHOS_ASSERT(thickness_);
#endif
  const auto & control_volumes = *(mesh_->control_volumes());
  auto c_data = control_volumes.getData();

#ifndef NDEBUG
  // Make sure control volumes and state still match.
  TEUCHOS_ASSERT_EQUALITY(2*c_data.size(), x_data.size());
#endif

  const auto thickness_values = thickness_->get_v(params);
  auto t_data = thickness_values.getData();
#ifndef NDEBUG
  TEUCHOS_ASSERT(control_volumes.getMap()->isSameAs(*thickness_values.getMap()));
#endif

  if (param_name.compare("g") == 0) {
    for (int k = 0; k < c_data.size(); k++) {
      // This assumes that "g" is not a parameter in either of the
      // potentials.
      double alpha = c_data[k] * t_data[k] *
        (x_data[2*k]*x_data[2*k] + x_data[2*k+1]*x_data[2*k+1]);
      // real and imaginary part
      f_data[2*k]   += alpha * x_data[2*k];
      f_data[2*k+1] += alpha * x_data[2*k+1];
    }
  } else {
    const auto scalar_potential_values =
      scalar_potential_->get_dvdp(params, param_name);
    auto s_data = scalar_potential_values.getData();
#ifndef NDEBUG
    TEUCHOS_ASSERT(control_volumes.getMap()->isSameAs(
          *scalar_potential_values.getMap()
          )
        );
#endif

    for (int k = 0; k < c_data.size(); k++) {
      const double alpha = c_data[k] * t_data[k] * s_data[k];
      // real and imaginary part
      f_data[2*k]   += alpha * x_data[2*k];
      f_data[2*k+1] += alpha * x_data[2*k+1];
    }
  }

  return;
}
// =============================================================================
double
nls::
inner_product(
    const Thyra::VectorBase<double> &phi,
    const Thyra::VectorBase<double> &psi
    ) const
{
  const auto & control_volumes = *mesh_->control_volumes();

  auto c_data = control_volumes.getData();

  size_t num_my_points = c_data.size();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*num_my_points, phi.space()->dim());
  TEUCHOS_ASSERT_EQUALITY(2*num_my_points, psi.space()->dim());
#endif

  double res = 0.0;
  for (size_t k = 0; k < num_my_points; k++) {
    double re_phi = Thyra::get_ele(phi, 2*k);
    double im_phi = Thyra::get_ele(phi, 2*k+1);
    double re_psi = Thyra::get_ele(psi, 2*k);
    double im_psi = Thyra::get_ele(psi, 2*k+1);
    res += c_data[k] * (re_phi*re_psi + im_phi*im_psi);
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
  //return globalRes / mesh_->control_volumes()->norm1();

  return 0.0;
}
// =============================================================================
double
nls::
gibbs_energy(const Thyra::VectorBase<double> &psi) const
{
  const auto & control_volumes = *mesh_->control_volumes();

  auto c_data = control_volumes.getData();

  size_t num_my_points = c_data.size();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*num_my_points, psi.space()->dim());
#endif

  double my_energy = 0.0;
  for (size_t k = 0; k < num_my_points; k++) {
    double re = Thyra::get_ele(psi, 2*k);
    double im = Thyra::get_ele(psi, 2*k+1);
    double alpha = re*re + im*im;
    my_energy -= c_data[k] * alpha * alpha;
  }

  // TODO
  //// Sum over all processors.
  //double globalEnergy;
  //TEUCHOS_ASSERT_EQUALITY(0, mesh_->getComm().SumAll(&my_energy, &globalEnergy, 1));

  //// normalize and return
  //return globalEnergy / mesh_->control_volumes()->norm1();

  return 0.0;
}
// =============================================================================
}  // namespace model_evaluator
}  // namespace nosh
