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
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_DefaultPreconditioner.hpp>
#include <Thyra_ModelEvaluatorBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_get_Epetra_Operator.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>

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
  const std::shared_ptr<Nosh::ParameterMatrix::Virtual> &keo,
  const std::shared_ptr<Nosh::ParameterMatrix::Virtual> &dKeoDP,
  const std::shared_ptr<const Nosh::ScalarField::Virtual> &scalarPotential,
  const double g,
  const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
  const std::shared_ptr<const Epetra_Vector> &initialX
) :
  mesh_(mesh),
  scalarPotential_(scalarPotential),
  thickness_(thickness),
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
  p_names_(Teuchos::null),
  nominalValues_(this->createInArgs())
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
  p_map_ = Teuchos::rcp(new Epetra_LocalMap(numParams, 0, mesh_->getComm()));
  Teuchos::RCP<Thyra::VectorBase<double>> p_init =
    Thyra::createMember(this->get_p_space(0));
  p_names_ = Teuchos::rcp(new Teuchos::Array<std::string>(numParams));
  int k = 0;
  for (auto it = params.begin(); it != params.end(); ++it) {
    (*p_names_)[k] = it->first;
    Thyra::set_ele(k, it->second, p_init());
    k++;
  }

  // set nominal values
  Teuchos::RCP<const Thyra::VectorBase<double>> xxx = Thyra::create_Vector(
        Teuchos::rcp(initialX), this->get_x_space()
        );
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
get_x_space() const
{
  // It is a bit of an assumption that x_ actually has this map, but as
  // Epetra_Vector::Map() only returns an Epetra_BlockMap which cannot be cast
  // into an Epetra_Map, this workaround is needed.
#ifndef NDEBUG
  TEUCHOS_ASSERT(mesh_);
#endif
  return Thyra::create_VectorSpace(Teuchos::rcp(
        mesh_->getComplexNonOverlapMap()
        ));
}
// ============================================================================
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
Nls::
get_f_space() const
{
#ifndef NDEBUG
  TEUCHOS_ASSERT(mesh_);
#endif
  return Thyra::create_VectorSpace(Teuchos::rcp(
        mesh_->getComplexNonOverlapMap()
        ));
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
  return Thyra::create_VectorSpace(p_map_);
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
  Teuchos::RCP<Epetra_Operator> jac =
      Teuchos::rcp(
        new Nosh::JacobianOperator(
          mesh_,
          scalarPotential_,
          thickness_,
          keo_
          )
        );

  return Thyra::nonconstEpetraLinearOp(jac, "Jacobian");
}
// =============================================================================
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
Nls::
get_W_factory() const
{
  Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "Belos");
  Teuchos::ParameterList & belosList =
    p->sublist("Linear Solver Types")
    .sublist("Belos");
  //belosList.set("Solver Type", "MINRES");
  belosList.set("Solver Type", "Pseudo Block GMRES");

  Teuchos::ParameterList & gmresList =
    belosList.sublist("Solver Types")
    .sublist("Pseudo Block GMRES");
  gmresList.set("Output Frequency", 10);
  gmresList.set("Output Style", 1);
  gmresList.set("Verbosity", 33);

  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double>>
    lowsFactory = builder.createLinearSolveStrategy("");

  lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  //// add preconditioner
  //Teuchos::RCP<Thyra::PreconditionerFactoryBase<double>>
  //  precFactory = Teuchos::rcp(new IfpackPreconditionerFactory());
  //if (precPL)
  //  precFactory->setParameterList(rcp(precPL,false));
  //lowsFactory->setPreconditionerFactory(precFactory,"Ifpack");

  return lowsFactory;
}
// =============================================================================
//Teuchos::RCP<Thyra::LinearOpWithSolveBase<double>>
//Nls::
//create_W() const
//{
//  Teuchos::RCP<const Thyra::LinearOpBase<double>> jac =
//    Thyra::epetraLinearOp(
//      Teuchos::rcp(
//        new Nosh::JacobianOperator(
//          mesh_,
//          scalarPotential_,
//          thickness_,
//          keo_
//          ),
//        "Jacobian"
//        ));
//
//   Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
//   // TODO set some solver parameters here
//   Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double>> lowsFactory =
//     linearSolverBuilder.createLinearSolveStrategy("");
//
//   return Thyra::linearOpWithSolve(*lowsFactory, jac);
//}
// =============================================================================
Teuchos::RCP<Thyra::PreconditionerBase<double>>
Nls::
create_W_prec() const
{
  const Teuchos::RCP<Epetra_Operator> keoPrec = Teuchos::rcp(
      new Nosh::KeoRegularized(
        mesh_,
        thickness_,
        keo_
        )
      );
  Teuchos::RCP<Thyra::LinearOpBase<double>> keoT =
    Thyra::nonconstEpetraLinearOp(keoPrec, "KEO");
  return Thyra::nonconstUnspecifiedPrec(keoT);
  //// bool is answer to: "Prec is already inverted?"
  //// This needs to be set to TRUE to make sure that the constructor of
  ////    NOX::Epetra::LinearSystemStratimikos
  //// chooses a user-defined preconditioner.
  //// Effectively, this boolean serves pretty well as a quirky switch for the
  //// preconditioner if Piro is used.
  //return Teuchos::rcp(
  //    new EpetraExt::ModelEvaluator::Preconditioner(keoPrec, true)
  //    );
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

  const Teuchos::RCP<const Thyra::VectorBase<double>> &x_in = inArgs.get_x();
#ifndef NDEBUG
  TEUCHOS_ASSERT(!x_in.is_null());
#endif
  // create corresponding epetra vector
  Teuchos::RCP<const Epetra_Comm> comm =
    Teuchos::rcpFromRef(mesh_->getComm());
  //Teuchos::RCP<const Epetra_Map> xmap =
  //  Thyra::get_Epetra_Map(*x_in->space(), comm);
  Teuchos::RCP<const Epetra_Vector> x_in_epetra =
    Thyra::get_Epetra_Vector(*mesh_->getComplexNonOverlapMap(), x_in);

  // Dissect inArgs.get_p(0) into parameter sublists.
  const Teuchos::RCP<const Thyra::VectorBase<double>> &p_in = inArgs.get_p(0);
#ifndef NDEBUG
  TEUCHOS_ASSERT(!p_in.is_null());
#endif

  Teuchos::RCP<const Epetra_Map> pmap =
    Thyra::get_Epetra_Map(*p_in->space(), comm);
  Teuchos::RCP<const Epetra_Vector> p_in_epetra =
    Thyra::get_Epetra_Vector(*pmap, p_in);

#ifndef NDEBUG
  // Make sure the parameters aren't NaNs.
  TEUCHOS_ASSERT(!p_in_epetra.is_null());
  for (int k = 0; k < p_in_epetra->MyLength(); k++) {
    TEUCHOS_ASSERT(!std::isnan((*p_in_epetra)[k]));
  }
#endif

  // Fill the parameters into a std::map.
  const Teuchos::RCP<const Teuchos::Array<std::string> > paramNames =
    this->get_p_names(0);
  std::map<std::string, double> params;
  for (int k = 0; k < p_in_epetra->MyLength(); k++) {
    params[(*paramNames)[k]] = (*p_in_epetra)[k];
    //std::cout << (*paramNames)[k] << " " << (*p_in)[k] << std::endl;
  }

  // compute F
  const Thyra::RCP<Thyra::VectorBase<double>> &f_out = outArgs.get_f();
  if (!f_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm1(*computeFTime_);
#endif

    //Teuchos::RCP<const Epetra_Map> fmap =
    //  Thyra::get_Epetra_Map(*f_out->space(), comm);
    Teuchos::RCP<Epetra_Vector> f_out_epetra =
      Thyra::get_Epetra_Vector(*mesh_->getComplexNonOverlapMap(), f_out);
    this->computeF_(
        *x_in_epetra,
        params,
        *f_out_epetra
        );
  }

  // Compute df/dp.
  const Thyra::ModelEvaluatorBase::DerivativeMultiVector<double> &derivMv =
    outArgs.get_DfDp(0).getDerivativeMultiVector();
  const Teuchos::RCP<Thyra::MultiVectorBase<double>> &dfdp_out =
    derivMv.getMultiVector();
  if (!dfdp_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm2(*computedFdpTime_);
#endif
    //Teuchos::RCP<const Epetra_Map> dfdpmap =
    //  Thyra::get_Epetra_Map(*dfdp_out->col(0)->space(), comm);
    Teuchos::RCP<Epetra_MultiVector> dfdp_out_epetra =
      Thyra::get_Epetra_MultiVector(
          *mesh_->getComplexNonOverlapMap(),
          dfdp_out
          );

    const int numAllParams = this->get_p_space(0)->dim();
    TEUCHOS_ASSERT_EQUALITY(
        numAllParams,
        dfdp_out_epetra->NumVectors()
        );
    // Compute all derivatives.
    for (int k = 0; k < numAllParams; k++) {
      this->computeDFDP_(
          *x_in_epetra,
          params,
          (*paramNames)[k],
          *(*dfdp_out_epetra)(k)
          );
    }
  }

  // Fill Jacobian.
  const Teuchos::RCP<Thyra::LinearOpBase<double>> & W_out = outArgs.get_W_op();
  if(!W_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm3(*fillJacobianTime_);
#endif
    Teuchos::RCP<Epetra_Operator> W_outE = Thyra::get_Epetra_Operator(*W_out);
    const Teuchos::RCP<Nosh::JacobianOperator> & jac =
      Teuchos::rcp_dynamic_cast<Nosh::JacobianOperator>(W_outE, true);
    jac->rebuild(params, *x_in_epetra);
  }

  // Fill preconditioner.
  const Teuchos::RCP<Thyra::PreconditionerBase<double>> & WPrec_out =
    outArgs.get_W_prec();
  if(!WPrec_out.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm4(*fillPreconditionerTime_);
#endif
    Teuchos::RCP<Epetra_Operator> WPrec_outE =
      Thyra::get_Epetra_Operator(
          *WPrec_out->getNonconstUnspecifiedPrecOp()
          );
    const Teuchos::RCP<Nosh::KeoRegularized> & keoPrec =
      Teuchos::rcp_dynamic_cast<Nosh::KeoRegularized>(WPrec_outE, true);
    keoPrec->rebuild(
        params,
        *x_in_epetra
        );
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
  TEUCHOS_ASSERT(mesh_);
  TEUCHOS_ASSERT(scalarPotential_);
  TEUCHOS_ASSERT(thickness_);
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
  TEUCHOS_ASSERT(mesh_);
  TEUCHOS_ASSERT(thickness_);
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
innerProduct(
    const Thyra::VectorBase<double> &phi,
    const Thyra::VectorBase<double> &psi
    ) const
{
  const Epetra_Vector &controlVolumes = *mesh_->getControlVolumes();

  int numMyPoints = controlVolumes.Map().NumMyPoints();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, phi.space()->dim());
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, psi.space()->dim());
#endif

  double res = 0.0;
  for (int k = 0; k < numMyPoints; k++) {
    double rePhi = Thyra::get_ele(phi, 2*k);
    double imPhi = Thyra::get_ele(phi, 2*k+1);
    double rePsi = Thyra::get_ele(psi, 2*k);
    double imPsi = Thyra::get_ele(psi, 2*k+1);
    res += controlVolumes[k] * (rePhi*rePsi + imPhi*imPsi);
  }

  // Sum over all processors.
  double globalRes;
  TEUCHOS_ASSERT_EQUALITY(0, mesh_->getComm().SumAll(&res, &globalRes, 1));

  // normalize and return
  return globalRes / mesh_->getDomainVolume();
}
// =============================================================================
double
Nls::
gibbsEnergy(const Thyra::VectorBase<double> &psi) const
{
  const Epetra_Vector &controlVolumes = *mesh_->getControlVolumes();

  int numMyPoints = controlVolumes.Map().NumMyPoints();
#ifndef NDEBUG
  TEUCHOS_ASSERT_EQUALITY(2*numMyPoints, psi.space()->dim());
#endif

  double myEnergy = 0.0;
  for (int k = 0; k < numMyPoints; k++) {
    double re = Thyra::get_ele(psi, 2*k);
    double im = Thyra::get_ele(psi, 2*k+1);
    double alpha = re*re + im*im;
    myEnergy -= controlVolumes[k] * alpha * alpha;
  }

  // Sum over all processors.
  double globalEnergy;
  TEUCHOS_ASSERT_EQUALITY(0, mesh_->getComm().SumAll(&myEnergy, &globalEnergy, 1));

  // normalize and return
  return globalEnergy / mesh_->getDomainVolume();
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
