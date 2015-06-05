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
Virtual(const std::shared_ptr<const Nosh::StkMesh> &mesh):
  Tpetra::CrsMatrix<double,int,int>(mesh->buildComplexGraph()),
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
setParameters(const std::map<std::string, double> &params)
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
    for (auto const &buildParam: buildParameters_) {
      if (buildParam.second != params.at(buildParam.first)) {
        needsRefill = true;
        break;
      }
    }
  }

  if (needsRefill) {
    this->refill_(params);
    for (auto &buildParam: buildParameters_) {
      buildParam.second = params.at(buildParam.first);
    }
  }

  return;
}
// ============================================================================
} // namespace ParameterMatrix
} // namespace Nosh
