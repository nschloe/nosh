// @HEADER
//
//    Builder class for the Laplace operator.
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
#ifndef NOSH_MATRIXBUILDER_LAPLACE_H
#define NOSH_MATRIXBUILDER_LAPLACE_H

#include <map>
#include <string>
#include <tuple>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Tpetra_CrsMatrix.hpp>

#include "Mesh.hpp"
#include "ParameterObject.hpp"

// forward declarations
namespace Nosh
{
class Mesh;
namespace ScalarField
{
class Virtual;
}
} // namespace Nosh

namespace Nosh
{
namespace ParameterMatrix
{
class Laplace: public Nosh::ParameterObject, public Tpetra::CrsMatrix<double,int,int>
{
public:
  Laplace(
      const std::shared_ptr<const Nosh::Mesh> &mesh,
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
  const std::shared_ptr<const Nosh::Mesh> mesh_;
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
