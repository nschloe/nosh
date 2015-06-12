// @HEADER
//
//    Nosh virtual model evaluator.
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
#ifndef NOSH_MODELEVALUATOR_VIRTUAL_H
#define NOSH_MODELEVALUATOR_VIRTUAL_H

// forward declarations
namespace Nosh
{
class Mesh;
}

#include <Thyra_ModelEvaluatorDefaultBase.hpp>

namespace Nosh
{
namespace ModelEvaluator
{
class Virtual : public Thyra::ModelEvaluatorDefaultBase<double>
{
public:
  Virtual();

  // Destructor
  virtual
  ~Virtual();

  virtual
  double
  innerProduct(
      const Thyra::VectorBase<double> &phi,
      const Thyra::VectorBase<double> &psi
      ) const = 0;

  virtual
  double
  norm(const Thyra::VectorBase<double> &psi) const;

  virtual
  double
  gibbsEnergy(const Thyra::VectorBase<double> &psi) const = 0;

  virtual
  const std::shared_ptr<const Nosh::Mesh>
  getMesh() const = 0;

protected:
private:
};
} // namespace ModelEvaluator
} // namespace Nosh

#endif // NOSH_MODELEVALUATOR_VIRTUAL_H
