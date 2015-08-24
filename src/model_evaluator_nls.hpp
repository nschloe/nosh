// @HEADER
//
//    Nosh model evaluator.
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
#ifndef NOSH_MODELEVALUATOR_NLS_H
#define NOSH_MODELEVALUATOR_NLS_H

// includes
#include <map>
#include <string>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_ParameterList.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "model_evaluator_base.hpp"

// forward declarations
namespace nosh
{
  class mesh;
  namespace scalar_field
  {
    class base;
  }
  namespace vector_field
  {
    class base;
  }
  namespace parameter_matrix
  {
    class keo;
    class DkeoDP;
  }
} // namespace nosh

namespace nosh
{
namespace model_evaluator
{
class nls : public base
{
public:
  nls (
    const std::shared_ptr<const nosh::mesh> &mesh,
    const std::shared_ptr<nosh::vector_field::base> &mvp,
    const std::shared_ptr<const nosh::scalar_field::base> &scalar_potential,
    const double g,
    const std::shared_ptr<const nosh::scalar_field::base> &thickness,
    const std::shared_ptr<const Tpetra::Vector<double,int,int>> &initial_x,
    const std::string & deriv_parameter
    );

  virtual
  ~nls();

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_x_space() const;

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_f_space() const;

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_p_space(int l) const;

  virtual
  Teuchos::RCP<const Teuchos::Array<std::string> >
  get_p_names(int l) const;

  virtual
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  get_g_space(int l) const;

  virtual
  Teuchos::ArrayView<const std::string>
  get_g_names(int j) const ;

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  getNominalValues() const
  {
    return nominal_values_;
  }

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  getLowerBounds() const;

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  getUpperBounds() const;

  //virtual
  //Teuchos::RCP<Thyra::LinearOpWithSolveBase<double>>
  //create_W() const;

  virtual
  Teuchos::RCP<Thyra::LinearOpBase<double>>
  create_W_op() const;

  virtual
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double>>
  get_W_factory() const;

  virtual
  Teuchos::RCP<Thyra::PreconditionerBase<double>>
  create_W_prec() const;

  virtual
  void
  reportFinalPoint(
      const Thyra::ModelEvaluatorBase::InArgs<double> &finalPoint,
      const bool wasSolved
      );

  virtual
  Thyra::ModelEvaluatorBase::InArgs<double>
  createInArgs() const;

public:
  double
  inner_product(
      const Thyra::VectorBase<double> &phi,
      const Thyra::VectorBase<double> &psi
      ) const;

  double
  gibbs_energy(const Thyra::VectorBase<double> &psi) const;

  const std::shared_ptr<const nosh::mesh>
  mesh() const
  {
    return mesh_;
  }

protected:

  virtual
  Thyra::ModelEvaluatorBase::OutArgs<double>
  createOutArgsImpl() const;

  virtual
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &in_args,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &out_args
      ) const;

private:
  Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
  createAlteredSpace() const;

  void
  compute_f_(
      const Tpetra::Vector<double,int,int> &x,
      const std::map<std::string, double> & params,
      Tpetra::Vector<double,int,int> &f_vec
      ) const;

  void
  computeDFDP_(
      const Tpetra::Vector<double,int,int> &x,
      const std::map<std::string, double> & params,
      const std::string & param_name,
      Tpetra::Vector<double,int,int> &f_vec
      ) const;

private:
  const std::shared_ptr<const nosh::mesh> mesh_;

  const std::shared_ptr<nosh::vector_field::base> mvp_;
  const std::shared_ptr<const nosh::scalar_field::base> scalar_potential_;
  const std::shared_ptr<const nosh::scalar_field::base> thickness_;

  const std::shared_ptr<nosh::parameter_matrix::keo> keo_;
  const std::shared_ptr<nosh::parameter_matrix::DkeoDP> dkeo_dp_;

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> eval_model_time_;
  const Teuchos::RCP<Teuchos::Time> compute_f_time_;
  const Teuchos::RCP<Teuchos::Time> compute_dfdp_time_;
  const Teuchos::RCP<Teuchos::Time> fill_jacobian_time_;
  const Teuchos::RCP<Teuchos::Time> fill_preconditioner_time_;
#endif

  Teuchos::RCP<Teuchos::FancyOStream> out_;

  Teuchos::RCP<const Tpetra::Map<int,int>> p_map_;
  Teuchos::RCP<Teuchos::Array<std::string> > p_names_;

  Thyra::ModelEvaluatorBase::InArgs<double> nominal_values_;

  const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> space_;
};
} // namespace model_evaluator
} // namespace nosh

#endif // NOSH_MODELEVALUATOR_NLS_H
