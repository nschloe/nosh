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

#ifndef NOSH_MATRIXBUILDER_KEO_H
#define NOSH_MATRIXBUILDER_KEO_H

#include <map>
#include <string>

#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif

#include <Tpetra_CrsMatrix.hpp>

#include <Eigen/Dense>

#include "nosh/StkMesh.hpp"
#include "nosh/ParameterObject.hpp"

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
} // namespace Nosh

namespace Nosh
{
namespace ParameterMatrix
{

class Keo: public Nosh::ParameterObject, public Tpetra::CrsMatrix<double,int,int>
{
public:
  Keo(
      const std::shared_ptr<const Nosh::StkMesh> &mesh,
      const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
      const std::shared_ptr<Nosh::VectorField::Virtual> &mvp
     );

  // Destructor.
  ~Keo();

  //! Gets the initial parameters from this module.
  virtual
  const std::map<std::string, double>
  getParameters() const;

protected:
private:
  void
  refill_(const std::map<std::string, double> & params);

  void
  buildAlphaCache_(
      const std::vector<edge> & edges,
      const std::vector<double> & edgeCoefficients
      ) const;

  double
  integrate1d_(
      const Nosh::VectorField::Virtual & f,
      const Eigen::Vector3d & x0,
      const Eigen::Vector3d & x1
      ) const;

private:
  const std::shared_ptr<const Nosh::StkMesh> mesh_;
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const std::shared_ptr<Teuchos::Time> keoFillTime_;
#endif
  const std::shared_ptr<const Nosh::ScalarField::Virtual> thickness_;
  const std::shared_ptr<Nosh::VectorField::Virtual> mvp_;

  mutable std::vector<double> alphaCache_;
  mutable bool alphaCacheUpToDate_;
};
} // namespace ParameterMatrix
} // namespace Nosh

#endif // NOSH_MATRIXBUILDER_KEO_H
