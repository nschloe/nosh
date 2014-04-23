// @HEADER
//
//    Nosh bordered model evaluator.
//    Copyright (C) 2012  Nico Schl\"omer
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

#include "nosh/ModelEvaluator_Bordered.hpp"

#include <string>
#include <vector>

#include "nosh/BorderingHelpers.hpp"
#include "nosh/BorderedOperator.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>

namespace Nosh
{
namespace ModelEvaluator
{
// ============================================================================
Bordered::
Bordered(const Teuchos::RCP<const Nosh::ModelEvaluator::Virtual> & modelEvaluator,
         const Teuchos::RCP<const Epetra_Vector> & initialBordering,
         const double lambdaInit
        ):
  innerModelEval_( modelEvaluator ),
  initialBordering_( initialBordering ),
  lambdaInit_( lambdaInit )
{
}
// ============================================================================
Bordered::
~Bordered()
{
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Bordered::
get_x_map() const
{
  return Nosh::BorderingHelpers::extendMapBy1(*innerModelEval_->get_x_map());
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Bordered::
get_f_map() const
{
  return Nosh::BorderingHelpers::extendMapBy1(*innerModelEval_->get_f_map());
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Bordered::
get_x_init() const
{
  const Teuchos::RCP<const Epetra_Vector> & inner_x_init =
    innerModelEval_->get_x_init();
  // Embed x_init into a vector of larger size.
  Teuchos::RCP<Epetra_Vector> out =
    Teuchos::rcp(new Epetra_Vector(*Nosh::BorderingHelpers::extendMapBy1(inner_x_init->Map())));

  Nosh::BorderingHelpers::merge(*inner_x_init,
                                &lambdaInit_,
                                *out);
  return out;
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Bordered::
get_p_init(int l) const
{
  return innerModelEval_->get_p_init(l);
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
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
Teuchos::RCP<Epetra_Operator>
Bordered::
create_W() const
{
  return Teuchos::rcp(new Nosh::BorderedOperator(innerModelEval_->create_W(),
                      Teuchos::rcp(new Epetra_Vector(*initialBordering_)),
                      Teuchos::rcp(new Epetra_Vector(*initialBordering_)),
                      0.0));
}
// =============================================================================
Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Bordered::
create_WPrec() const
{
  // Extract the inner preconditioner and fill a Bordered Operator with it.
  Teuchos::RCP<Epetra_Operator> borderedPrec =
    Teuchos::rcp(new Nosh::BorderedOperator(innerModelEval_->create_WPrec()->PrecOp,
                 Teuchos::rcp(new Epetra_Vector(*initialBordering_)),
                 Teuchos::rcp(new Epetra_Vector(*initialBordering_)),
                 0.0));

  // bool is answer to: "Prec is already inverted?"
  // This needs to be set to TRUE to make sure that the constructor of
  //    NOX::Epetra::LinearSystemStratimikos
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
evalModel(const InArgs &inArgs,
          const OutArgs &outArgs
         ) const
{
  // First, dissect x_in into vector and bordering.
  const Teuchos::RCP<const Epetra_Vector> &x_in = inArgs.get_x();
#ifndef NDEBUG
  TEUCHOS_ASSERT(!x_in.is_null());
#endif
  const Teuchos::RCP<Epetra_Vector> inner_x_in =
    Teuchos::rcp(new Epetra_Vector(*innerModelEval_->get_x_map()));
  double lambda[1];
  Nosh::BorderingHelpers::dissect(*x_in, *inner_x_in, lambda);

  // Get i*x. This assumes a particular data layout in x_in.
  const Teuchos::RCP<Epetra_Vector> ix =
    Teuchos::rcp(new Epetra_Vector(inner_x_in->Map()));
  for (int k=0; k<ix->Map().NumMyElements()/2; k++) {
    (*ix)[2*k] = - (*x_in)[2*k+1];
    (*ix)[2*k+1] = (*x_in)[2*k];
  }

  // Copy over the args for use in innerModelEval.
  InArgs inner_inArgs = inArgs;
  inner_inArgs.set_x(inner_x_in);

  OutArgs inner_outArgs = outArgs;

  const Teuchos::RCP<const Epetra_Vector> & bordering =
    //initialBordering_;
    ix;

  // Compute F(x).
  const Teuchos::RCP<Epetra_Vector> &f_out = outArgs.get_f();
  if ( !f_out.is_null() ) {
    // Create new temporary f_out.
    const Teuchos::RCP<Epetra_Vector> inner_f_out =
      Teuchos::rcp(new Epetra_Vector(*innerModelEval_->get_f_map()));

    inner_outArgs.set_f(inner_f_out);
    innerModelEval_->evalModel(inner_inArgs, inner_outArgs);
    // Add lambda * x0.
    TEUCHOS_ASSERT_EQUALITY(0, inner_f_out->Update(lambda[0], *bordering, 1.0));
    // Append <psi0, x> to f_out.
    double r[1];
    TEUCHOS_ASSERT_EQUALITY(0, bordering->Dot(*inner_x_in, r));
    //r = lambda;
    Nosh::BorderingHelpers::merge(*inner_f_out, r, *f_out);
  }

  // Compute df/dp.
  const EpetraExt::ModelEvaluator::DerivativeMultiVector &derivMv =
    outArgs.get_DfDp(0).getDerivativeMultiVector();
  const Teuchos::RCP<Epetra_MultiVector> &dfdp_out =
    derivMv.getMultiVector();
  if ( !dfdp_out.is_null() ) {
    // Create temporary DerivativeMultiVector inner_dfdp_out.
    const int numParams = derivMv.getParamIndexes().length();
    const Teuchos::RCP<Epetra_MultiVector> inner_dfdp_out =
      Teuchos::rcp(new Epetra_MultiVector(*innerModelEval_->get_f_map(),
                                          numParams));
    const EpetraExt::ModelEvaluator::DerivativeMultiVector innerDerivMv(inner_dfdp_out,
        derivMv.getOrientation(),
        derivMv.getParamIndexes());
    inner_outArgs.set_DfDp(0, innerDerivMv);
    innerModelEval_->evalModel(inner_inArgs, inner_outArgs);
    // Append last entry and merge into dfdp_out.
    std::vector<double> r(numParams);
    for (int k=0; k<numParams; k++)
      r[k] = 0.0;
    Nosh::BorderingHelpers::merge(*inner_dfdp_out, &r[0], *dfdp_out);
  }

  // Fill Jacobian.
  const Teuchos::RCP<Epetra_Operator> & W_out = outArgs.get_W();
  if( !W_out.is_null() ) {
    const Teuchos::RCP<Nosh::BorderedOperator> & borderedW =
      Teuchos::rcp_dynamic_cast<Nosh::BorderedOperator>(W_out, true);

    // Fill inner Jacobian.
    inner_outArgs.set_W(borderedW->getInnerOperator());
    innerModelEval_->evalModel(inner_inArgs, inner_outArgs);

    // Reset bordering.
    borderedW->resetBordering(bordering, bordering, 0.0);
  }

  // Fill preconditioner.
  const Teuchos::RCP<Epetra_Operator> & WPrec_out = outArgs.get_WPrec();
  if( !WPrec_out.is_null() ) {
    const Teuchos::RCP<Nosh::BorderedOperator> & borderedPrec =
      Teuchos::rcp_dynamic_cast<Nosh::BorderedOperator>(WPrec_out, true);

    // Fill inner preconditioner.
    inner_outArgs.set_WPrec(borderedPrec->getInnerOperator());
    innerModelEval_->evalModel(inner_inArgs, inner_outArgs);

    // Reset bordering.
    borderedPrec->resetBordering(bordering, bordering, 0.0);
  }

  return;
}
// =============================================================================
double
Bordered::
innerProduct(const Epetra_Vector &phi,
             const Epetra_Vector &psi
            ) const
{
  return innerModelEval_->innerProduct(phi, psi);
}
// =============================================================================
double
Bordered::
gibbsEnergy(const Epetra_Vector &psi) const
{
  return innerModelEval_->gibbsEnergy(psi);
}
// =============================================================================
const Teuchos::RCP<const Nosh::StkMesh>
Bordered::
getMesh() const
{
  return innerModelEval_->getMesh();
}
// =============================================================================
} // namespace ModelEvaluator
} // namespace Nosh
