// @HEADER
//
//    Regularized kinetic energy operator.
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

#ifndef NOSH_KEOPRECONDITIONER_H
#define NOSH_KEOPRECONDITIONER_H

#include <map>
#include <string>

#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Teuchos_FancyOStream.hpp>

#include <MueLu_EpetraOperator.hpp>

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

namespace Nosh
{
class KeoRegularized : public Epetra_Operator
{
public:
  KeoRegularized(
      const std::shared_ptr<const Nosh::StkMesh> &mesh,
      const std::shared_ptr<const Nosh::ScalarField::Virtual> &thickness,
      const std::shared_ptr<const Nosh::ParameterMatrix::Virtual> &matrix
      );

  // Destructor.
  ~KeoRegularized();

  virtual int
  SetUseTranspose(bool UseTranspose);

  virtual int
  Apply(const Epetra_MultiVector &X,
        Epetra_MultiVector &Y
       ) const;

  virtual int
  ApplyInverse(const Epetra_MultiVector &X,
               Epetra_MultiVector &Y
              ) const;

  virtual double
  NormInf() const;

  virtual const char *
  Label() const;

  virtual bool
  UseTranspose() const;

  virtual bool
  HasNormInf() const;

  virtual const Epetra_Comm &
  Comm() const;

  virtual const Epetra_Map &OperatorDomainMap() const;

  virtual const Epetra_Map &OperatorRangeMap() const;

public:
  void
  rebuild(const std::map<std::string, double> & params,
          const Epetra_Vector &psi
         );

protected:
private:
  const std::shared_ptr<const Epetra_Vector>
  getAbsPsiSquared_(const Epetra_Vector &psi);

  void
  rebuildInverse_();

private:
  bool useTranspose_;

  const std::shared_ptr<const Nosh::StkMesh> mesh_;
  const std::shared_ptr<const Nosh::ScalarField::Virtual> thickness_;

  const std::shared_ptr<Nosh::ParameterMatrix::Virtual> regularizedKeo_;

  Teuchos::RCP<MueLu::EpetraOperator> MueluPrec_;

  const int numCycles_;

#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const std::shared_ptr<Teuchos::Time> timerRebuild0_;
  const std::shared_ptr<Teuchos::Time> timerRebuild1_;
#endif

  Teuchos::RCP<Teuchos::FancyOStream> out_;
};
} // namespace Nosh

#endif // NOSH_KEOPRECONDITIONER_H
