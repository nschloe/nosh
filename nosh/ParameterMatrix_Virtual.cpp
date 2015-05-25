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

#include "nosh/ParameterMatrix_Virtual.hpp"

#include "nosh/StkMesh.hpp"

namespace Nosh
{
namespace ParameterMatrix
{
// ============================================================================
Virtual::
Virtual(const Teuchos::RCP<const Nosh::StkMesh> &mesh):
  Epetra_FECrsMatrix(Copy, mesh->buildComplexGraph()),
  mesh_(mesh),
  buildParameters_()
{
}
// ============================================================================
Virtual::
~Virtual()
{
}
// ============================================================================
void
Virtual::
refill(const std::map<std::string, double> &params)
{
  // Cache the construction of the matrix.
  // This is useful because in the continuation context, the matrix is called a
  // number of times with the same arguments (in computeF, getJacobian(), and
  // getPreconditioner().
  bool needsRefill;
  if (buildParameters_.empty()) {
    needsRefill = true;
  } else {
    needsRefill = false;
    for (auto it = buildParameters_.begin();
         it != buildParameters_.end();
         ++it) {
      // Check if it->first is in params at all and if their values are equal.
      std::map<std::string, double>::const_iterator it2 =
        params.find(it->first);
      TEUCHOS_ASSERT(it2 != params.end());
      if (it2->second != it->second) {
        needsRefill = true;
        break;
      }
    }
  }

  if (needsRefill) {
    this->refill_(params);
    // Reset build parameters.
    for (auto it = buildParameters_.begin();
         it != buildParameters_.end();
         ++it) {
      std::map<std::string, double>::const_iterator it2 =
        params.find(it->first);
      TEUCHOS_ASSERT(it2 != params.end());
      it->second = it2->second;
    }
  }

  return;
}
// ============================================================================
} // namespace ParameterMatrix
} // namespace Nosh
