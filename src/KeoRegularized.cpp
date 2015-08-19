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

#include "KeoRegularized.hpp"

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

#include "ScalarField_Virtual.hpp"
#include "VectorField_Virtual.hpp"
#include "ParameterMatrix_Keo.hpp"
#include "Mesh.hpp"

// =============================================================================
// some typdefs for Belos
using ST = double;
using MV = Tpetra::MultiVector<double,int,int>;
using OP = Tpetra::Operator<double,int,int>;
using MVT = Belos::MultiVecTraits<ST, MV>;
using OPT = Belos::OperatorTraits<ST, MV, OP>;
// =============================================================================
namespace Nosh
{
// =============================================================================
KeoRegularized::
KeoRegularized(
    const std::shared_ptr<const Nosh::Mesh> &mesh,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
    const std::shared_ptr<Nosh::VectorField::Virtual> &mvp
    ):
  mesh_(mesh),
  thickness_(thickness),
  // It wouldn't strictly be necessary to initialize regularizedKeo_ with the
  // proper graph here as matrixBuilder_'s cache will override the matrix later
  // on anyways. Keep it, though, as it doesn't waste any memory and is in the
  // spirit of the Trilinos::ModelEvaluator which asks for allocation of memory
  // at one point and filling it with meaningful values later on.
  regularizedKeo_(
      std::make_shared<Nosh::ParameterMatrix::Keo>(mesh, thickness, mvp)
      ),
  MueluPrec_(Teuchos::null),
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  timerRebuild0_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: KeoRegularized::rebuild::MueLu init"
        )),
  timerRebuild1_(Teuchos::TimeMonitor::getNewTimer(
        "Nosh: KeoRegularized::rebuild::MueLu rebuild"
        )),
#endif
  numCycles_(1)
{
}
// =============================================================================
KeoRegularized::
~KeoRegularized()
{
}
// =============================================================================
//int
//KeoRegularized::
//applyInverse(
//    const Tpetra::MultiVector<double,int,int> &X,
//    Tpetra::MultiVector<double,int,int> &Y
//    ) const
//{
//#ifndef NDEBUG
//  TEUCHOS_ASSERT(regularizedKeo_->DomainMap().SameAs(X.getMap()));
//  TEUCHOS_ASSERT(regularizedKeo_->RangeMap().SameAs(Y.getMap()));
//#endif
//  // (K +  g * 2*|psi|) * X
//  return regularizedKeo_->Apply(X, Y);
//}
// =============================================================================
void
KeoRegularized::
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

  if (numCycles_ == 1) {
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
    belosList.set("Maximum Iterations", numCycles_);
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
        Teuchos::rcp(regularizedKeo_),
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
KeoRegularized::
getDomainMap() const
{
  return regularizedKeo_->getDomainMap();
}
// =============================================================================
Teuchos::RCP<const Tpetra::Map<int,int>>
KeoRegularized::
getRangeMap() const
{
  return regularizedKeo_->getRangeMap();
}
// =============================================================================
void
KeoRegularized::
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
  // worry too much about this until memory contrains get tight.
  regularizedKeo_->setParameters(params);

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
    //TEUCHOS_ASSERT(regularizedKeo_.RowMap().SameAs(absPsiSquared->Map()));
//#endif
    //Tpetra::Vector<double,int,int> diag(regularizedKeo_.RowMap());
    //TEUCHOS_ASSERT_EQUALITY(0, regularizedKeo_.ExtractDiagonalCopy(diag));
    //TEUCHOS_ASSERT_EQUALITY(0, diag.Update(g*2.0, *absPsiSquared, 1.0));
    //TEUCHOS_ASSERT_EQUALITY(0, regularizedKeo_.ReplaceDiagonalValues(diag));
    //
    const auto & controlVolumes = *(mesh_->getControlVolumes());
    const auto thicknessValues = thickness_->getV(params);
#ifndef NDEBUG
    TEUCHOS_ASSERT(controlVolumes.getMap()->isSameAs(*thicknessValues.getMap()));
#endif
    auto cData = controlVolumes.getData();
    auto tData = thicknessValues.getData();
    auto xData = x.getData();
#ifndef NDEBUG
    TEUCHOS_ASSERT_EQUALITY(cData.size(), tData.size());
    TEUCHOS_ASSERT_EQUALITY(2*tData.size(), xData.size());
#endif
    Teuchos::Tuple<int,2> idx;
    Teuchos::Tuple<double,2> vals;
    regularizedKeo_->resumeFill();
    for (int k = 0; k < cData.size(); k++) {
      const double alpha = g * cData[k] * tData[k]
        * 2.0 * (xData[2*k]*xData[2*k] + xData[2*k+1]*xData[2*k+1]);

      const double beta = g * cData[k] * tData[k]
                          * (2.0 * xData[2*k] * xData[2*k+1]);

      const double gamma = g * cData[k] * tData[k]
        * (xData[2*k]*xData[2*k] - xData[2*k+1]*xData[2*k+1]);

      // TODO check if the indices are correct here
      idx[0] = 2*k;
      idx[1] = 2*k + 1;
      vals[0] = alpha + gamma;
      vals[1] = beta;
      TEUCHOS_ASSERT_EQUALITY(
          2,
          regularizedKeo_->sumIntoLocalValues(2*k, idx, vals)
          );
      vals[0] = beta;
      vals[1] = alpha - gamma;
      TEUCHOS_ASSERT_EQUALITY(
          2,
          regularizedKeo_->sumIntoLocalValues(2*k+1, idx, vals)
          );
    }
  }
  regularizedKeo_->fillComplete();

  this->rebuildInverse_();
  return;
}
// =============================================================================
//const std::shared_ptr<const Tpetra::Vector<double,int,int>>
//KeoRegularized::
//getAbsPsiSquared_(const Tpetra::Vector<double,int,int> &psi)
//{
//#ifndef NDEBUG
//  TEUCHOS_ASSERT(!mesh_.is_null());
//  TEUCHOS_ASSERT(!thickness_.is_null());
//#endif
//  const std::shared_ptr<Tpetra::Vector<double,int,int>> absPsiSquared =
//    Teuchos::rcp(new Tpetra::Vector<double,int,int>(psi.getMap()));
//
//  const std::shared_ptr<const Tpetra::Vector<double,int,int>> &controlVolumes =
//    mesh_->getControlVolumes();
//  int numMyPoints = controlVolumes->MyLength();
//  for (int k=0; k<numMyPoints; k++)
//  {
//    double alpha = (*controlVolumes)[k] * thickness_->getV(k)
//                 * (psi[2*k]*psi[2*k] + psi[2*k+1]*psi[2*k+1]);
//    (*absPsiSquared)[2*k] = alpha;
//    (*absPsiSquared)[2*k+1] = alpha;
//  }
//  return absPsiSquared;
//}
// =============================================================================
void
KeoRegularized::
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
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int>> rKeoRcp =
      Teuchos::rcp(regularizedKeo_);

    MueluPrec_ = MueLu::CreateTpetraPreconditioner(
        rKeoRcp,
        params
        );
  } else {
#ifdef NOSH_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor tm(*timerRebuild1_);
#endif
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int>> rKeoRcp =
      Teuchos::rcp(regularizedKeo_);
    MueLu::ReuseTpetraPreconditioner(
        rKeoRcp,
        *MueluPrec_
        );
  }

  return;
}
// =============================================================================
}  // namespace Nosh
