// @HEADER
//
//    Nosh Jacobian operator.
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

#ifndef NOSH_JACOBIANOPERATOR_H
#define NOSH_JACOBIANOPERATOR_H

#include <map>
#include <string>

#include <Tpetra_Vector.hpp>
#include <Tpetra_Operator.hpp>
#include <Teuchos_RCP.hpp>

// forward declarations
namespace Nosh
{
  class Mesh;
  namespace ParameterMatrix
  {
    class Keo;
  }
  namespace ScalarField
  {
    class Virtual;
  }
} // namespace Nosh

namespace Nosh
{

class JacobianOperator : public Tpetra::Operator<double,int,int>
{
public:
  JacobianOperator(
      const std::shared_ptr<const Nosh::Mesh> &mesh,
      const std::shared_ptr<const Nosh::ScalarField::Virtual> &scalarPotential,
      const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
      const std::shared_ptr<Nosh::ParameterMatrix::Keo> &keo
      );

  // Destructor.
  ~JacobianOperator ();

  virtual void
  apply(
      const Tpetra::MultiVector<double,int,int> &X,
      Tpetra::MultiVector<double,int,int> &Y,
      Teuchos::ETransp mode,
      double alpha,
      double beta
      ) const;

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>> getDomainMap() const;

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>> getRangeMap() const;

public:
  void
  rebuild(
      const std::map<std::string, double> params,
      const Tpetra::Vector<double,int,int> & current_X
      );

protected:
private:
  void
  rebuildDiags_(
      const std::map<std::string, double> params,
      const Tpetra::Vector<double,int,int> &current_X
      );

private:
  const std::shared_ptr<const Nosh::Mesh> mesh_;
  const std::shared_ptr<const Nosh::ScalarField::Virtual> scalarPotential_;
  const std::shared_ptr<const Nosh::ScalarField::Virtual> thickness_;

  const std::shared_ptr<Nosh::ParameterMatrix::Keo> keo_;
  Tpetra::Vector<double,int,int> diag0_;
  Tpetra::Vector<double,int,int> diag1b_;
};
} // namespace Nosh

#endif // NOSH_JACOBIANOPERATOR_H
