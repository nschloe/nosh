// @HEADER
//
//    Nosh virtual model evaluator.
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
#ifndef NOSH_MODELEVALUATOR_BASE_H
#define NOSH_MODELEVALUATOR_BASE_H

// forward declarations
namespace nosh
{
class mesh;
}

#include <Thyra_ModelEvaluatorDefaultBase.hpp>

namespace nosh
{
namespace model_evaluator
{
class base : public Thyra::ModelEvaluatorDefaultBase<double>
{
public:
  base()
  {
  }

  // Destructor
  virtual
  ~base()
  {
  }

  virtual
  double
  inner_product(
      const Thyra::VectorBase<double> &phi,
      const Thyra::VectorBase<double> &psi
      ) const = 0;

  virtual
  double
  norm(const Thyra::VectorBase<double> &psi) const
  {
    return sqrt(this->inner_product(psi, psi));
  }

  virtual
  double
  gibbs_energy(const Thyra::VectorBase<double> &psi) const = 0;

  virtual
  const std::shared_ptr<const nosh::mesh>
  mesh() const = 0;

protected:
private:
};
} // namespace model_evaluator
} // namespace nosh

#endif // NOSH_MODELEVALUATOR_BASE_H
