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
#include <tuple>

#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>
#ifdef NOSH_TEUCHOS_TIME_MONITOR
#include <Teuchos_Time.hpp>
#endif
#include <Teuchos_Array.hpp>
#include <Epetra_FECrsGraph.h>

#include <stk_mesh/base/Entity.hpp>

#include "nosh/MatrixBuilder_Virtual.hpp"

typedef std::tuple<stk::mesh::Entity, stk::mesh::Entity> edge;

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
namespace MatrixBuilder
{

class Keo: public Virtual
{
public:
  Keo(const Teuchos::RCP<const Nosh::StkMesh> &mesh,
      const Teuchos::RCP<const Nosh::ScalarField::Virtual> &thickness,
      const Teuchos::RCP<const Nosh::VectorField::Virtual> &mvp
    );

  // Destructor.
  ~Keo();

  virtual
  void
  apply(
      const std::map<std::string, double> & params,
      const Epetra_Vector &X,
      Epetra_Vector &Y
      ) const;

  virtual
  void
  applyDKDp(
      const std::map<std::string, double> & params,
      const std::string & paramName,
      const Epetra_Vector &X,
      Epetra_Vector &Y
      ) const;

  virtual
  void
  fill(
      Epetra_FECrsMatrix &matrix,
      const std::map<std::string, double> & params
      ) const;

//! Gets the initial parameters from this module.
  virtual
  const std::map<std::string, double>
  getInitialParameters() const;

protected:
private:
  void
  fillKeo_(
      Epetra_FECrsMatrix &keoMatrix,
      const std::map<std::string, double> & params,
      void (Keo::*filler)(const int, const std::map<std::string, double>&, double*) const
      ) const;

  void
  fillerRegular_(
      const int k,
      const std::map<std::string, double> & params,
      double * v
      ) const;

  void
  fillerDp_(const int k,
            const std::map<std::string, double> & params,
            double * v
            ) const;

  void
  buildAlphaCache_(const Teuchos::Array<edge> & edges,
                   const Teuchos::ArrayRCP<const double> &edgeCoefficients
                   ) const;

private:
#ifdef NOSH_TEUCHOS_TIME_MONITOR
  const Teuchos::RCP<Teuchos::Time> keoFillTime_;
  const Teuchos::RCP<Teuchos::Time> buildKeoGraphTime_;
#endif
  const Teuchos::RCP<const Nosh::ScalarField::Virtual> thickness_;
  const Teuchos::RCP<const Nosh::VectorField::Virtual> mvp_;

  mutable Epetra_FECrsMatrix keoCache_;
  mutable std::map<std::string, double> keoBuildParameters_;
  mutable Epetra_FECrsMatrix keoDpCache_;


  mutable Teuchos::ArrayRCP<double> alphaCache_;
  mutable bool alphaCacheUpToDate_;
  mutable std::string paramName_;
};
} // namespace MatrixBuilder
} // namespace Nosh

#endif // NOSH_MATRIXBUILDER_KEO_H
