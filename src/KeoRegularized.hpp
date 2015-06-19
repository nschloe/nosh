// @HEADER
//
//    Regularized kinetic energy operator.
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

#ifndef NOSH_KEOPRECONDITIONER_H
#define NOSH_KEOPRECONDITIONER_H

#include <map>
#include <string>

#include <Tpetra_Vector.hpp>
#include <Tpetra_Operator.hpp>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Teuchos_FancyOStream.hpp>

#include <MueLu_TpetraOperator.hpp>

// forward declarations
namespace Nosh
{
  class Mesh;
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
    class Keo;
  }
} // namespace Nosh

namespace Nosh
{
class KeoRegularized : public Tpetra::Operator<double,int,int>
{
public:
  KeoRegularized(
      const std::shared_ptr<const Nosh::Mesh> &mesh,
      const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
      const std::shared_ptr<Nosh::VectorField::Virtual> &mvp
      );

  // Destructor.
  ~KeoRegularized();

  virtual void
  apply(
      const Tpetra::MultiVector<double,int,int> &X,
      Tpetra::MultiVector<double,int,int> &Y,
      Teuchos::ETransp mode,
      double alpha,
      double beta
      ) const;

  //virtual int
  //applyInverse(const Tpetra::MultiVector<double,int,int> &X,
  //             Tpetra::MultiVector<double,int,int> &Y
  //            ) const;

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>> getDomainMap() const;

  virtual
  Teuchos::RCP<const Tpetra::Map<int,int>> getRangeMap() const;

public:
  void
  rebuild(
      const std::map<std::string, double> & params,
      const Tpetra::Vector<double,int,int> &psi
      );

protected:
private:
  const std::shared_ptr<const Tpetra::Vector<double,int,int>>
  getAbsPsiSquared_(const Tpetra::Vector<double,int,int> &psi);

  void
  rebuildInverse_();

private:

  const std::shared_ptr<const Nosh::Mesh> mesh_;
  const std::shared_ptr<const Nosh::ScalarField::Virtual> thickness_;

  const std::shared_ptr<Nosh::ParameterMatrix::Keo> regularizedKeo_;

  Teuchos::RCP<MueLu::TpetraOperator<double,int,int>> MueluPrec_;

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> timerRebuild0_;
  const Teuchos::RCP<Teuchos::Time> timerRebuild1_;
#endif

  const int numCycles_;
};
} // namespace Nosh

#endif // NOSH_KEOPRECONDITIONER_H
