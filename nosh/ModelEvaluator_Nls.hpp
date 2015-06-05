// @HEADER
//
//    Nosh model evaluator.
//    Copyright (C) 2010--2012  Nico Schl\"omer
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

#include "nosh/ModelEvaluator_Virtual.hpp"

// forward declarations
namespace Nosh
{
  class StkMesh;
  namespace ScalarField
  {
    class Virtual;
  }
  namespace VectorField
  {
    class Virtual;
  }
  namespace ParameterMatrix
  {
    class Virtual;
  }
} // namespace Nosh

namespace Nosh
{
namespace ModelEvaluator
{
class Nls : public Virtual
{
public:
  //! Constructor without initial guess.
  Nls (
    const std::shared_ptr<const Nosh::StkMesh> &mesh,
    const std::shared_ptr<Nosh::VectorField::Virtual> &mvp,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &scalarPotential,
    const double g,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
    const std::shared_ptr<const Tpetra::Vector<double,int,int>> &initialX
    );

  // Destructor
  virtual
  ~Nls();

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
  Thyra::ModelEvaluatorBase::InArgs<double>
  getNominalValues() const;

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
  innerProduct(
      const Thyra::VectorBase<double> &phi,
      const Thyra::VectorBase<double> &psi
      ) const;

  double
  gibbsEnergy(const Thyra::VectorBase<double> &psi) const;

  const std::shared_ptr<const Nosh::StkMesh>
  getMesh() const;

protected:

  virtual
  Thyra::ModelEvaluatorBase::OutArgs<double>
  createOutArgsImpl() const;

  virtual
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs
      ) const;

private:
  void
  computeF_(
      const Tpetra::Vector<double,int,int> &x,
      const std::map<std::string, double> & params,
      Tpetra::Vector<double,int,int> &FVec
      ) const;

  void
  computeDFDP_(
      const Tpetra::Vector<double,int,int> &x,
      const std::map<std::string, double> & params,
      const std::string & paramName,
      Tpetra::Vector<double,int,int> &FVec
      ) const;

private:
  const std::shared_ptr<const Nosh::StkMesh> mesh_;

  const std::shared_ptr<Nosh::VectorField::Virtual> mvp_;
  const std::shared_ptr<const Nosh::ScalarField::Virtual> scalarPotential_;
  const std::shared_ptr<const Nosh::ScalarField::Virtual> thickness_;

  const std::shared_ptr<Nosh::ParameterMatrix::Virtual> keo_;
  const std::shared_ptr<Nosh::ParameterMatrix::Virtual> dKeoDP_;

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const std::shared_ptr<Teuchos::Time> evalModelTime_;
  const std::shared_ptr<Teuchos::Time> computeFTime_;
  const std::shared_ptr<Teuchos::Time> computedFdpTime_;
  const std::shared_ptr<Teuchos::Time> fillJacobianTime_;
  const std::shared_ptr<Teuchos::Time> fillPreconditionerTime_;
#endif

  Teuchos::RCP<Teuchos::FancyOStream> out_;

  Teuchos::RCP<const Tpetra::Map<int,int>> p_map_;
  Teuchos::RCP<Teuchos::Array<std::string> > p_names_;

  Thyra::ModelEvaluatorBase::InArgs<double> nominalValues_;
};
} // namespace ModelEvaluator
} // namespace Nosh

#endif // NOSH_MODELEVALUATOR_NLS_H
