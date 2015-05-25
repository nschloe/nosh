// @HEADER
//
//    Builder class for the kinetic energy operator.
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

#ifndef NOSH_MATRIXBUILDER_DKEODP_H
#define NOSH_MATRIXBUILDER_DKEODP_H

#include <map>
#include <string>
#include <tuple>

#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Teuchos_Array.hpp>
#include <Epetra_FECrsGraph.h>

#include <stk_mesh/base/Entity.hpp>

#include "nosh/ParameterMatrix_Virtual.hpp"

#include "nosh/StkMesh.hpp"

// forward declarations
namespace Nosh
{
namespace ScalarField
{
class Virtual;
}
namespace VectorField
{
class Virtual;
}
} // namespace Nosh

namespace Nosh
{
namespace ParameterMatrix
{

class DKeoDP: public Virtual
{
public:
  DKeoDP(
      const Teuchos::RCP<const Nosh::StkMesh> &mesh,
      const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness,
      const Teuchos::RCP<const Nosh::VectorField::Virtual> &mvp,
      const std::string & paramName
      );

  // Destructor.
  ~DKeoDP();

  virtual
  Teuchos::RCP<Virtual>
  clone() const;

  //! Gets the initial parameters from this module.
  virtual
  const std::map<std::string, double>
  getInitialParameters() const;

protected:
private:
  void
  refill_(const std::map<std::string, double> & params);

  void
  buildAlphaCache_(
      const Teuchos::Array<edge> & edges,
      const Teuchos::ArrayRCP<const double> &edgeCoefficients
      ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> keoFillTime_;
#endif
  const Teuchos::RCP<const Nosh::ScalarField::Virtual> thickness_;
  const Teuchos::RCP<const Nosh::VectorField::Virtual> mvp_;

  mutable Teuchos::ArrayRCP<double> alphaCache_;
  mutable bool alphaCacheUpToDate_;
  const std::string paramName_;
};
} // namespace ParameterMatrix
} // namespace Nosh

#endif // NOSH_MATRIXBUILDER_DKEODP_H
