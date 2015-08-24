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
#ifndef NOSH_MODELEVALUATOR_BORDERED_H
#define NOSH_MODELEVALUATOR_BORDERED_H
// -----------------------------------------------------------------------------
// includes
#include <string>

#include <Teuchos_RCP.hpp>

#include "model_evaluator_base.hpp"
// -----------------------------------------------------------------------------
namespace nosh
{
namespace model_evaluator
{
class Bordered : public nosh::model_evaluator::base
{
public:
  //! Constructor without initial guess.
  Bordered (
      const std::shared_ptr<const nosh::model_evaluator::base> & model_eval,
      const std::shared_ptr<const Tpetra::Vector<double,int,int>> & initialBordering,
      const double lambdaInit
      );

  // Destructor
  virtual
  ~Bordered();

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>>
  get_x_map() const;

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>>
  get_f_map() const;

  virtual
  Teuchos::RCP<const Tpetra::Vector<double,int,int>>
  get_x_init() const;

  virtual
  Teuchos::RCP<const Tpetra::Vector<double,int,int>>
  get_p_init(int l) const;

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>>
  get_p_map(int l) const;

  virtual
  Teuchos::RCP<const Teuchos::Array<std::string> >
  get_p_names(int l) const;

  virtual
  Teuchos::RCP<Tpetra::Operator<double,int,int>>
  create_W() const;

  virtual
  Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
  create_WPrec() const;

  virtual
  InArgs
  createInArgs() const;

  virtual
  OutArgs
  createOutArgs() const;

  virtual
  void
  eval_mdel(
      const InArgs &in_args,
      const OutArgs &out_args
      ) const;

public:
  virtual
  double
  inner_product(const Tpetra::Vector<double,int,int> &phi,
               const Tpetra::Vector<double,int,int> &psi
             ) const;

  virtual
  double
  gibbs_energy(const Tpetra::Vector<double,int,int> &psi) const;

  virtual
  const std::shared_ptr<const nosh::mesh>
  mesh() const;

protected:
private:
  const std::shared_ptr<const nosh::model_evaluator::base> innerModelEval_;
  const std::shared_ptr<const Tpetra::Vector<double,int,int>> initialBordering_;
  const double lambdaInit_;
};
} // namespace model_evaluator
} // namespace nosh

#endif // NOSH_MODELEVALUATOR_BORDERED_H
