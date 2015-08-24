// @HEADER
//
//    Regularized kinetic energy operator.
//    Copyright (C) 2010--2012  Nico Schlömer
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

#include "keo_regularized.hpp"

#include <map>
#include <string>

#include <Teuchos_Comm.hpp>
#include <Tpetra_Vector.hpp>
#include <Epetra_FECrsMatrix.h>

#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosOperatorT.hpp>
#include <BelosMueLuAdapter.hpp>

#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

#include "scalar_field_base.hpp"
#include "vector_field_base.hpp"
#include "parameter_matrix_keo.hpp"
#include "mesh.hpp"

// =============================================================================
// some typdefs for Belos
using ST = double;
using MV = Tpetra::MultiVector<double,int,int>;
using OP = Tpetra::Operator<double,int,int>;
using MVT = Belos::MultiVecTraits<ST, MV>;
using OPT = Belos::OperatorTraits<ST, MV, OP>;
// =============================================================================
namespace nosh
{
// =============================================================================
keo_regularized::
keo_regularized(
    const std::shared_ptr<const nosh::mesh> &mesh,
    const std::shared_ptr<const nosh::scalar_field::base> &thickness,
    const std::shared_ptr<nosh::vector_field::base> &mvp
    ):
  mesh_(mesh),
  thickness_(thickness),
  // It wouldn't strictly be necessary to initialize regularizedkeo_ with the
  // proper graph here as matrixBuilder_'s cache will override the matrix later
  // on anyways. Keep it, though, as it doesn't waste any memory and is in the
  // spirit of the Trilinos::ModelEvaluator which asks for allocation of memory
  // at one point and filling it with meaningful values later on.
  regularizedkeo_(
      std::make_shared<nosh::parameter_matrix::keo>(mesh, thickness, mvp)
      ),
  MueluPrec_(Teuchos::null),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  timerRebuild0_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: keo_regularized::rebuild::MueLu init"
        )),
  timerRebuild1_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: keo_regularized::rebuild::MueLu rebuild"
        )),
#endif
  num_cycles_(1)
{
}
// =============================================================================
keo_regularized::
~keo_regularized()
{
}
// =============================================================================
//int
//keo_regularized::
//applyInverse(
//    const Tpetra::MultiVector<double,int,int> &X,
//    Tpetra::MultiVector<double,int,int> &Y
//    ) const
//{
//#ifndef NDEBUG
//  TEUCHOS_ASSERT(regularizedkeo_->DomainMap().SameAs(X.getMap()));
//  TEUCHOS_ASSERT(regularizedkeo_->RangeMap().SameAs(Y.getMap()));
//#endif
//  // (K +  g * 2*|psi|) * X
//  return regularizedkeo_->Apply(X, Y);
//}
// =============================================================================
void
keo_regularized::
apply(
    const Tpetra::MultiVector<double,int,int> &X,
    Tpetra::MultiVector<double,int,int> &Y,
    Teuchos::ETransp mode,
    double alpha,
    double beta
    ) const
{
  TEUCHOS_ASSERT_EQUALITY(mode, Teuchos::NO_TRANS);
  TEUCHOS_ASSERT_EQUALITY(alpha, 1.0);
  TEUCHOS_ASSERT_EQUALITY(beta, 0.0);

#ifndef NDEBUG
  TEUCHOS_ASSERT(!MueluPrec_.is_null());
#endif

  if (num_cycles_ == 1) {
    // Just apply one (inverse) AMG cycle.
    return MueluPrec_->apply(X, Y);
  } else {
    // Belos part
    Teuchos::ParameterList belosList;
    // Relative convergence tolerance requested.
    // Set this to 0 and adapt the maximum number of iterations. This way, the
    // preconditioner always does exactly the same thing (namely maxIter PCG
    // iterations) and is independent of X. This avoids mathematical
    // difficulties.
    belosList.set("Convergence Tolerance", 0.0);
    belosList.set("Maximum Iterations", num_cycles_);
//     belosList.set("Verbosity",
//                   Belos::Errors +
//                   Belos::Warnings +
//                   Belos::TimingDetails +
//                   Belos::StatusTestDetails
//                 );
//     belosList.set("Output Frequency", 10);
    belosList.set("Verbosity", Belos::Errors + Belos::Warnings);

    // Make sure to have a solid initial guess.
    // Belos, for example, does not initialize Y before passing it here.
    Y.putScalar(0.0);

    // Construct an unpreconditioned linear problem instance.
    auto Xptr = Teuchos::rcpFromRef(X);
    auto Yptr = Teuchos::rcpFromRef(Y);
    Belos::LinearProblem<double, MV, OP> problem(
        Teuchos::rcp(regularizedkeo_),
        Yptr,
        Xptr
        );
    // Make sure the problem sets up correctly.
    TEUCHOS_ASSERT(problem.setProblem());

    // add preconditioner
    // TODO recheck
    // Create the Belos preconditioned operator from the preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    //Teuchos::RCP<Belos::EpetraPrecOp> mueluPrec =
    //  Teuchos::rcp(new Belos::EpetraPrecOp(MueluPrec_));
    problem.setLeftPrec(MueluPrec_);

    // Create an iterative solver manager.
    Belos::PseudoBlockCGSolMgr<double, MV, OP> newSolver(
        Teuchos::rcp(&problem, false),
        Teuchos::rcp(&belosList, false)
        );

    // Perform "solve".
    Belos::ReturnType ret = newSolver.solve();

    TEUCHOS_ASSERT_EQUALITY(ret, Belos::Converged);

    return;
  }
}
// =============================================================================
Teuchos::RCP<const Tpetra::Map<int,int>>
keo_regularized::
getDomainMap() const
{
  return regularizedkeo_->getDomainMap();
}
// =============================================================================
Teuchos::RCP<const Tpetra::Map<int,int>>
keo_regularized::
getRangeMap() const
{
  return regularizedkeo_->getRangeMap();
}
// =============================================================================
void
keo_regularized::
rebuild(
    const std::map<std::string, double> & params,
    const Tpetra::Vector<double,int,int> & x
    )
{
  // Copy over the matrix.
  // This is necessary as we don't apply AMG to K, but to K + g*2|psi|^2.
  // A possible work-around this copy would be to define the matrixBuilder's
  // matrix as K+2*g*|psi|^2 in the first place, and make sure that 2*g*|psi|^2
  // is taken away again wherever needed (e.g., the Jacobian).  This would
  // introduce the additional complication of having KEO depend on psi, and
  // would likely lead to some confusion in the rest of the code.  Hence, don't
  // worry too much about this until memory constrains get tight.
  regularizedkeo_->set_parameters(params);

  const double g = params.at("g");
  // Add 2*g*|psi|^2 to the diagonal.
  if (g > 0.0) {
    // To be added to the diagonal blocks:
    //
    // [alpha + gamma, beta         ]
    // [beta,          alpha - gamma].
    //
    // We could also ahead and only add alpha to the diagonal, i.e.,
    //
    //const std::shared_ptr<const Tpetra::Vector<double,int,int>> absPsiSquared =
    //  this->getAbsPsiSquared_(x);
//#ifndef NDEBUG
    //TEUCHOS_ASSERT(regularizedkeo_.RowMap().SameAs(absPsiSquared->Map()));
//#endif
    //Tpetra::Vector<double,int,int> diag(regularizedkeo_.RowMap());
    //TEUCHOS_ASSERT_EQUALITY(0, regularizedkeo_.ExtractDiagonalCopy(diag));
    //TEUCHOS_ASSERT_EQUALITY(0, diag.Update(g*2.0, *absPsiSquared, 1.0));
    //TEUCHOS_ASSERT_EQUALITY(0, regularizedkeo_.ReplaceDiagonal_values(diag));
    //
    const auto & control_volumes = *(mesh_->control_volumes());
    const auto thickness_values = thickness_->get_v(params);
#ifndef NDEBUG
    TEUCHOS_ASSERT(control_volumes.getMap()->isSameAs(*thickness_values.getMap()));
#endif
    auto c_data = control_volumes.getData();
    auto t_data = thickness_values.getData();
    auto x_data = x.getData();
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(c_data.size(), t_data.size());
    TEUCHOS_ASSERT_EQUALITY(2*t_data.size(), x_data.size());
#endif
    Teuchos::Tuple<int,2> idx;
    Teuchos::Tuple<double,2> vals;
    regularizedkeo_->resumeFill();
    for (int k = 0; k < c_data.size(); k++) {
      const double alpha = g * c_data[k] * t_data[k]
        * 2.0 * (x_data[2*k]*x_data[2*k] + x_data[2*k+1]*x_data[2*k+1]);

      const double beta = g * c_data[k] * t_data[k]
                          * (2.0 * x_data[2*k] * x_data[2*k+1]);

      const double gamma = g * c_data[k] * t_data[k]
        * (x_data[2*k]*x_data[2*k] - x_data[2*k+1]*x_data[2*k+1]);

      // TODO check if the indices are correct here
      idx[0] = 2*k;
      idx[1] = 2*k + 1;
      vals[0] = alpha + gamma;
      vals[1] = beta;
      TEUCHOS_ASSERT_EQUALITY(
          2,
          regularizedkeo_->sumIntoLocalValues(2*k, idx, vals)
          );
      vals[0] = beta;
      vals[1] = alpha - gamma;
      TEUCHOS_ASSERT_EQUALITY(
          2,
          regularizedkeo_->sumIntoLocalValues(2*k+1, idx, vals)
          );
    }
  }
  regularizedkeo_->fillComplete();

  this->rebuildInverse_();
  return;
}
// =============================================================================
//const std::shared_ptr<const Tpetra::Vector<double,int,int>>
//keo_regularized::
//getAbsPsiSquared_(const Tpetra::Vector<double,int,int> &psi)
//{
//#ifndef NDEBUG
//  TEUCHOS_ASSERT(!mesh_.is_null());
//  TEUCHOS_ASSERT(!thickness_.is_null());
//#endif
//  const std::shared_ptr<Tpetra::Vector<double,int,int>> absPsiSquared =
//    Teuchos::rcp(new Tpetra::Vector<double,int,int>(psi.getMap()));
//
//  const std::shared_ptr<const Tpetra::Vector<double,int,int>> &control_volumes =
//    mesh_->control_volumes();
//  int num_my_points = control_volumes->MyLength();
//  for (int k=0; k<num_my_points; k++)
//  {
//    double alpha = (*control_volumes)[k] * thickness_->get_v(k)
//                 * (psi[2*k]*psi[2*k] + psi[2*k+1]*psi[2*k+1]);
//    (*absPsiSquared)[2*k] = alpha;
//    (*absPsiSquared)[2*k+1] = alpha;
//  }
//  return absPsiSquared;
//}
// =============================================================================
void
keo_regularized::
rebuildInverse_()
{
  if (MueluPrec_.is_null()) {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm(*timerRebuild0_);
#endif
    Teuchos::ParameterList params;
    params.set("number of equations", 2);
    params.set("reuse: type", "full");
    // See
    // <http://trilinos.org/wordpress/wp-content/uploads/2015/05/MueLu_Users_Guide_Trilinos12_0.pdf>
    // for recommended settings.
    //params.set("smoother: type", "Chebyshev");

    //params.set("max levels", 10);
    //params.set("increasing or decreasing", "increasing");
    //params.set("aggregation: type", "Uncoupled");
    //params.set("aggregation: threshold", 0.0);
    //params.set("smoother: sweeps", 3);
    //params.set("smoother: pre or post", "both");
    //params.set("coarse: type", "Amesos-KLU");

    // For some reason, we need to rcp explicitly. Otherwise, the call to
    // MueLu::CreateTpetraPreconditioner will complain about unmatching types.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int>> rkeoRcp =
      Teuchos::rcp(regularizedkeo_);

    MueluPrec_ = MueLu::CreateTpetraPreconditioner(
        rkeoRcp,
        params
        );
  } else {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm(*timerRebuild1_);
#endif
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int>> rkeoRcp =
      Teuchos::rcp(regularizedkeo_);
    MueLu::ReuseTpetraPreconditioner(
        rkeoRcp,
        *MueluPrec_
        );
  }

  return;
}
// =============================================================================
}  // namespace nosh
