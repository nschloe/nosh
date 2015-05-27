// @HEADER
//
//    Query routines for a vector potential with explicitly given values.
//    Copyright (C) 2011, 2012  Nico Schl\"omer
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
#ifndef NOSH_SCALARFIELD_EXPLICITVALUES_H_
#define NOSH_SCALARFIELD_EXPLICITVALUES_H_

#include <map>
#include <string>

#include <Teuchos_RCP.hpp>
#include <Epetra_MultiVector.h>

#include "nosh/ScalarField_Virtual.hpp"
#include "nosh/StkMesh.hpp"

namespace Nosh
{
namespace ScalarField
{
class ExplicitValues : public Virtual
{
public:
  ExplicitValues(
      const Nosh::StkMesh &mesh,
      const std::string &fieldName
      );

  virtual
  ~ExplicitValues();

//! Get parameter names and initial values.
  virtual
  const std::map<std::string, double>
  getParameters() const;

  virtual
  const Epetra_Vector
  getV(const std::map<std::string, double> & params) const;

  virtual
  const Epetra_Vector
  getdVdP(
      const std::map<std::string, double> & params,
      const std::string & paramName
      ) const;

protected:
private:
  const Teuchos::RCP<const Epetra_Vector> nodeValues_;
};
} // namespace ScalarField
} // namespace Nosh
#endif // NOSH_SCALARFIELD_EXPLICITVALUES_H_
