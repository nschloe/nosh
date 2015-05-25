// @HEADER
//
//    Virtual class for matrix constructors.
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
#ifndef NOSH_MATRIXBUILDER_VIRTUAL
#define NOSH_MATRIXBUILDER_VIRTUAL

#include <map>
#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Epetra_Comm.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>

#include <stk_mesh/base/Entity.hpp>

// forward declarations
namespace Nosh
{
class StkMesh;
}

namespace Nosh
{
namespace ParameterMatrix
{
class Virtual: public Epetra_FECrsMatrix
{
public:
  Virtual(const Teuchos::RCP<const Nosh::StkMesh> &mesh);

  // Destructor.
  virtual
  ~Virtual();

  //// https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Virtual_Constructor
  //virtual
  //Teuchos::RCP<Virtual>
  //create() const = 0;

  virtual
  Teuchos::RCP<Virtual>
  clone() const = 0;

  //! Fill the matrix with the parameter entries as given in params.
  //! Includes some caching logic for params.
  virtual
  void
  refill(const std::map<std::string, double> &params);

  //! Get parameter map with their initial values.
  virtual
  const std::map<std::string, double>
  getInitialParameters() const = 0;

protected:
  //! Fill the matrix with the parameter entries as given in params.
  virtual
  void
  refill_(const std::map<std::string, double> &params) = 0;

protected:
  const Teuchos::RCP<const Nosh::StkMesh> mesh_;

private:
  const Epetra_FECrsGraph
  buildGraph_(const Nosh::StkMesh &mesh);

private:
  std::map<std::string, double> buildParameters_;
};
} // namespace ParameterMatrix
} // namespace Nosh

#endif // NOSH_MATRIXBUILDER_VIRTUAL
