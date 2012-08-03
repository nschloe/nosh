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

#include "Nosh_ModelEvaluator_Bordered.hpp"

namespace Nosh {
namespace ModelEvaluator {
// ============================================================================
Bordered::
Bordered(const Teuchos::RCP<const Nosh::ModelEvaluator::Virtual> & modelEvaluator):
  innerModelEval_( modelEvaluator )
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
  return innerModelEval_->get_x_map();
}
// ============================================================================
Teuchos::RCP<const Epetra_Map>
Bordered::
get_f_map() const
{
  return innerModelEval_->get_f_map();
}
// ============================================================================
Teuchos::RCP<const Epetra_Vector>
Bordered::
get_x_init() const
{
  return innerModelEval_->get_x_init();
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
  return innerModelEval_->create_W();
}
// =============================================================================
Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Bordered::
create_WPrec() const
{
  return innerModelEval_->create_WPrec();
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
  innerModelEval_->evalModel(inArgs, outArgs);
  return;
}
// =============================================================================
Teuchos::RCP<const Epetra_Vector>
Bordered::
get_p_latest() const
{
  return innerModelEval_->get_p_latest();
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
normalizedScaledL2Norm(const Epetra_Vector &psi) const
{
  return innerModelEval_->normalizedScaledL2Norm(psi);
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
