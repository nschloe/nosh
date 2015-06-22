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
#ifndef NOSH_LAPLACE_H
#define NOSH_LAPLACE_H

#include <map>
#include <string>
#include <tuple>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include "LinearOperator.hpp"
#include "Mesh.hpp"
#include "ParameterObject.hpp"

// forward declarations
namespace Nosh
{
class Mesh;
} // namespace Nosh

namespace Nosh
{
class Laplace:
  public LinearOperator
{
public:
  Laplace(
      const std::shared_ptr<const Nosh::Mesh> & mesh,
      const std::shared_ptr<const Nosh::DirichletBoundaryConditions> & _bcs
      );

  // Destructor.
  ~Laplace();

protected:
private:
  void
  fill_();

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> fillTime_;
#endif
};
} // namespace Nosh

#endif // NOSH_LAPLACE_H
