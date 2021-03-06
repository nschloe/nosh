// @HEADER
//
//    Nosh bordered model evaluator.
//    Copyright (C) 2012  Nico Schlömer
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

#include "ModelEvaluator_Bordered.hpp"

#include <string>
#include <vector>

#include "BorderingHelpers.hpp"
#include "BorderedOperator.hpp"

#include <Teuchos::Comm<int>.h>
#include <Tpetra::Map<int,int>.h>
#include <Epetra_Import.h>

#include <Teuchos_RCPStdSharedPtrConversions.hpp>

namespace nosh
{
namespace model_evaluator
{
// ============================================================================
Bordered::
Bordered(const std::shared_ptr<const nosh::model_evaluator::base> & model_evaluator,
         const std::shared_ptr<const Tpetra::Vector<double,int,int>> & initialBordering,
         const double lambdaInit
       ):
  innerModelEval_(model_evaluator),
  initialBordering_(initialBordering),
  lambdaInit_(lambdaInit)
{
}
// ============================================================================
Bordered::
~Bordered()
{
}
// ============================================================================
Teuchos::RCP<const Tpetra::Map<int,int>>
Bordered::
get_x_map() const
{
  return Teuchos::rcp(
      nosh::BorderingHelpers::extendMapBy1(*innerModelEval_->get_x_map())
      );
}
// ============================================================================
Teuchos::RCP<const Tpetra::Map<int,int>>
Bordered::
get_f_map() const
{
  return Teuchos::rcp(
      nosh::BorderingHelpers::extendMapBy1(*innerModelEval_->get_f_map())
      );
}
// ============================================================================
Teuchos::RCP<const Tpetra::Vector<double,int,int>>
Bordered::
get_x_init() const
{
  const Teuchos::RCP<const Tpetra::Vector<double,int,int>> & inner_x_init =
    innerModelEval_->get_x_init();
  // Embed x_init into a vector of larger size.
  Teuchos::RCP<Tpetra::Vector<double,int,int>> out =
    Teuchos::rcp(new Tpetra::Vector<double,int,int>(*nosh::BorderingHelpers::extendMapBy1(inner_x_init->Map())));

  nosh::BorderingHelpers::merge(*inner_x_init,
                                &lambdaInit_,
                                *out);
  return out;
}
// ============================================================================
Teuchos::RCP<const Tpetra::Vector<double,int,int>>
Bordered::
get_p_init(int l) const
{
  return innerModelEval_->get_p_init(l);
}
// ============================================================================
Teuchos::RCP<const Tpetra::Map<int,int>>
Bordered::
get_p_map(int l) const
{
  return innerModelEval_->get_p_map(l);
}
// ============================================================================
Teuchos::RCP<const Teuchos::Array<std::string> >
Bordered::
get_p_names(int l) const
{
  return innerModelEval_->get_p_names(l);
}
// =============================================================================
Teuchos::RCP<Tpetra::Operator<double,int,int>>
Bordered::
create_W() const
{
  return Teuchos::rcp(
      new nosh::BorderedOperator(
        Teuchos::get_shared_ptr(innerModelEval_->create_W()),
        *initialBordering_,
        *initialBordering_,
        0.0
        )
      );
}
// =============================================================================
Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Bordered::
create_WPrec() const
{
  // Extract the inner preconditioner and fill a Bordered Operator with it.
  Teuchos::RCP<Tpetra::Operator<double,int,int>> borderedPrec =
    Teuchos::rcp(
        new nosh::BorderedOperator(
          Teuchos::get_shared_ptr(innerModelEval_->create_WPrec()->PrecOp),
          *initialBordering_,
          *initialBordering_,
          0.0
          )
        );

  // bool is answer to: "Prec is already inverted?"
  // This needs to be set to TRUE to make sure that the constructor of
  //    NOX::Epetra::linear_systemStratimikos
  // chooses a user-defined preconditioner.
  // Effectively, this boolean serves pretty well as a quirky switch for the
  // preconditioner if Piro is used.
  return Teuchos::rcp(new EpetraExt::ModelEvaluator::Preconditioner(borderedPrec,
                      true));
  //false));
}
// ============================================================================
EpetraExt::ModelEvaluator::InArgs
Bordered::
createInArgs() const
{
  return innerModelEval_->createInArgs();
}
// ============================================================================
EpetraExt::ModelEvaluator::OutArgs
Bordered::
createOutArgs() const
{
  return innerModelEval_->createOutArgs();
}
// ============================================================================
void
Bordered::
eval_mdel(const InArgs &in_args,
          const OutArgs &out_args
        ) const
{
  // First, dissect x_in into vector and bordering.
  const Teuchos::RCP<const Tpetra::Vector<double,int,int>> &x_in = in_args.get_x();
#ifndef NDEBUG
  TEUCHOS_ASSERT(!x_in.is_null());
#endif
  const Teuchos::RCP<Tpetra::Vector<double,int,int>> inner_x_in =
    Teuchos::rcp(new Tpetra::Vector<double,int,int>(*innerModelEval_->get_x_map()));
  double lambda[1];
  nosh::BorderingHelpers::dissect(*x_in, *inner_x_in, lambda);

  // Get i*x. This assumes a particular data layout in x_in.
  Tpetra::Vector<double,int,int> ix(inner_x_in->Map());
  for (int k = 0; k < ix.getMap().NumMyElements()/2; k++) {
    ix[2*k] = - (*x_in)[2*k+1];
    ix[2*k+1] = (*x_in)[2*k];
  }

  // Copy over the args for use in innerModelEval.
  InArgs inner_in_args = in_args;
  inner_in_args.set_x(inner_x_in);

  OutArgs inner_out_args = out_args;

  const Tpetra::Vector<double,int,int> & bordering = ix;

  // Compute F(x).
  const Teuchos::RCP<Tpetra::Vector<double,int,int>> &f_out = out_args.get_f();
  if (!f_out.is_null()) {
    // Create new temporary f_out.
    const Teuchos::RCP<Tpetra::Vector<double,int,int>> inner_f_out =
      Teuchos::rcp(new Tpetra::Vector<double,int,int>(*innerModelEval_->get_f_map()));

    inner_out_args.set_f(inner_f_out);
    innerModelEval_->eval_mdel(inner_in_args, inner_out_args);
    // Add lambda * x0.
    TEUCHOS_ASSERT_EQUALITY(0, inner_f_out->Update(lambda[0], bordering, 1.0));
    // Append <psi0, x> to f_out.
    double r[1];
    TEUCHOS_ASSERT_EQUALITY(0, bordering.Dot(*inner_x_in, r));
    //r = lambda;
    nosh::BorderingHelpers::merge(*inner_f_out, r, *f_out);
  }

  // Compute df/dp.
  const EpetraExt::ModelEvaluator::DerivativeMultiVector &derivMv =
    out_args.get_DfDp(0).getDerivativeMultiVector();
  const Teuchos::RCP<Tpetra::MultiVector<double,int,int>> &dfdp_out =
    derivMv.multi_vector();
  if (!dfdp_out.is_null()) {
    // Create temporary DerivativeMultiVector inner_dfdp_out.
    const int numParams = derivMv.get_paramIndexes().length();
    const Teuchos::RCP<Tpetra::MultiVector<double,int,int>> inner_dfdp_out =
      Teuchos::rcp(new Tpetra::MultiVector<double,int,int>(*innerModelEval_->get_f_map(),
                                          numParams));
    const EpetraExt::ModelEvaluator::DerivativeMultiVector innerDerivMv(inner_dfdp_out,
        derivMv.getOrientation(),
        derivMv.get_paramIndexes());
    inner_out_args.set_DfDp(0, innerDerivMv);
    innerModelEval_->eval_mdel(inner_in_args, inner_out_args);
    // Append last entry and merge into dfdp_out.
    std::vector<double> r(numParams);
    for (int k = 0; k < numParams; k++)
      r[k] = 0.0;
    nosh::BorderingHelpers::merge(*inner_dfdp_out, &r[0], *dfdp_out);
  }

  // Fill Jacobian.
  const Teuchos::RCP<Tpetra::Operator<double,int,int>> & W_out = out_args.get_W();
  if(!W_out.is_null()) {
    const Teuchos::RCP<nosh::BorderedOperator> & borderedW =
      Teuchos::rcp_dynamic_cast<nosh::BorderedOperator>(W_out, true);

    // Fill inner Jacobian.
    inner_out_args.set_W(Teuchos::rcp(borderedW->getInnerOperator()));
    innerModelEval_->eval_mdel(inner_in_args, inner_out_args);

    // Reset bordering.
    borderedW->resetBordering(bordering, bordering, 0.0);
  }

  // Fill preconditioner.
  const Teuchos::RCP<Tpetra::Operator<double,int,int>> & WPrec_out = out_args.get_WPrec();
  if(!WPrec_out.is_null()) {
    const Teuchos::RCP<nosh::BorderedOperator> & borderedPrec =
      Teuchos::rcp_dynamic_cast<nosh::BorderedOperator>(WPrec_out, true);

    // Fill inner preconditioner.
    inner_out_args.set_WPrec(Teuchos::rcp(borderedPrec->getInnerOperator()));
    innerModelEval_->eval_mdel(inner_in_args, inner_out_args);

    // Reset bordering.
    borderedPrec->resetBordering(bordering, bordering, 0.0);
  }

  return;
}
// =============================================================================
double
Bordered::
inner_product(const Tpetra::Vector<double,int,int> &phi,
             const Tpetra::Vector<double,int,int> &psi
           ) const
{
  return innerModelEval_->inner_product(phi, psi);
}
// =============================================================================
double
Bordered::
gibbs_energy(const Tpetra::Vector<double,int,int> &psi) const
{
  return innerModelEval_->gibbs_energy(psi);
}
// =============================================================================
const std::shared_ptr<const nosh::mesh>
Bordered::
mesh() const
{
  return innerModelEval_->mesh();
}
// =============================================================================
} // namespace model_evaluator
} // namespace nosh
