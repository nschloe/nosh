// @HEADER
//
//    Builder class for the Laplace operator.
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
#ifndef NOSH_MATRIXBUILDER_LAPLACE_H
#define NOSH_MATRIXBUILDER_LAPLACE_H

#include <map>
#include <string>
#include <tuple>

#include <Tpetra_Operator.hpp>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Tpetra_CrsGraph.hpp>

#include <stk_mesh/base/Entity.hpp>

#include "nosh/ParameterMatrix_Virtual.hpp"
#include "nosh/StkMesh.hpp"

// forward declarations
namespace Nosh
{
class StkMesh;
namespace ScalarField
{
class Virtual;
}
} // namespace Nosh

namespace Nosh
{
namespace ParameterMatrix
{
class Laplace: public Virtual
{
public:
  Laplace(
      const std::shared_ptr<const Nosh::StkMesh> &mesh,
      const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness
      );

  // Destructor.
  ~Laplace();

  //! Gets the parameter with their initial values.
  virtual
  const std::map<std::string, double>
  getParameters() const;

protected:
private:
  void
  fill_();

  void
  buildAlphaCache_(
      const std::vector<edge> & edges,
      const std::vector<double> & edgeCoefficients
      ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const std::shared_ptr<Teuchos::Time> fillTime_;
#endif
  const std::shared_ptr<const Nosh::ScalarField::Virtual> thickness_;

  mutable std::vector<double> alphaCache_;
  mutable bool alphaCacheUpToDate_;
};
} // namespace ParameterMatrix
} // namespace Nosh

#endif // NOSH_MATRIXBUILDER_LAPLACE_H
