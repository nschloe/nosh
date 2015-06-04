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
// -----------------------------------------------------------------------------
// includes
#include <map>
#include <string>

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "nosh/ModelEvaluator_Virtual.hpp"
// -----------------------------------------------------------------------------
// forward declarations
namespace Nosh
{
class StkMesh;
namespace ScalarField
{
class Virtual;
}
namespace ParameterMatrix
{
class Virtual;
}
} // namespace Nosh

class Epetra_CrsGraph;
class Epetra_LocalMap;
// -----------------------------------------------------------------------------
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
    const std::shared_ptr<Nosh::ParameterMatrix::Virtual> &keo,
    const std::shared_ptr<Nosh::ParameterMatrix::Virtual> &dKeoDP,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &scalarPotential,
    const double g,
    const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
    const std::shared_ptr<const Epetra_Vector> &initialX
    );

  // Destructor
  virtual
  ~Nls();

  virtual
  Teuchos::RCP<const Epetra_Map>
  get_x_map() const;

  virtual
  Teuchos::RCP<const Epetra_Map>
  get_f_map() const;

  virtual
  Teuchos::RCP<const Epetra_Vector>
  get_x_init() const;

  virtual
  Teuchos::RCP<const Epetra_Vector>
  get_p_init(int l) const;

  virtual
  Teuchos::RCP<const Epetra_Map>
  get_p_map(int l) const;

  virtual
  Teuchos::RCP<const Teuchos::Array<std::string> >
  get_p_names(int l) const;

  virtual
  Teuchos::RCP<Epetra_Operator>
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
  evalModel(
      const InArgs &inArgs,
      const OutArgs &outArgs
      ) const;

public:
  double
  innerProduct(
      const Epetra_Vector &phi,
      const Epetra_Vector &psi
      ) const;

  double
  gibbsEnergy(const Epetra_Vector &psi) const;

  const std::shared_ptr<const Nosh::StkMesh>
  getMesh() const;

private:
  void
  computeF_(const Epetra_Vector &x,
            const std::map<std::string, double> & params,
            Epetra_Vector &FVec
          ) const;

  void
  computeDFDP_(const Epetra_Vector &x,
               const std::map<std::string, double> & params,
               const std::string & paramName,
               Epetra_Vector &FVec
             ) const;

protected:
private:
  const std::shared_ptr<const Nosh::StkMesh> mesh_;

  const std::shared_ptr<const Nosh::ScalarField::Virtual> scalarPotential_;

  const std::shared_ptr<const Nosh::ScalarField::Virtual> thickness_;

  const Teuchos::RCP<const Epetra_Vector> x_init_;

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

  Teuchos::RCP<Epetra_LocalMap> p_map_;
  Teuchos::RCP<Epetra_Vector> p_init_;
  Teuchos::RCP<Teuchos::Array<std::string> > p_names_;
};
} // namespace ModelEvaluator
} // namespace Nosh

#endif // NOSH_MODELEVALUATOR_NLS_H
