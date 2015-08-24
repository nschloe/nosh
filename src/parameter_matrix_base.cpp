// @HEADER
//
//    base class for matrix constructors.
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

#include "parameter_matrix_base.hpp"

#include "mesh.hpp"

namespace nosh
{
namespace parameter_matrix
{
// ============================================================================
base::
base(const std::shared_ptr<const nosh::mesh> &mesh):
  Tpetra::CrsMatrix<double,int,int>(mesh->build_complex_graph()),
  mesh_(mesh),
  build_parameters_()
{
}
// ============================================================================
base::
~base()
{
}
// ============================================================================
void
base::
set_parameters(const std::map<std::string, double> &params)
{
  // Cache the construction of the matrix.
  // This is useful because in the continuation context, the matrix is called a
  // number of times with the same arguments (in compute_f, getJacobian(), and
  // get_preconditioner().
  bool needs_refill;
  if (build_parameters_.empty()) {
    needs_refill = true;
  } else {
    needs_refill = false;
    for (auto const &build_param: build_parameters_) {
      if (build_param.second != params.at(build_param.first)) {
        needs_refill = true;
        break;
      }
    }
  }

  if (needs_refill) {
    this->refill_(params);
    for (auto &build_param: build_parameters_) {
      build_param.second = params.at(build_param.first);
    }
  }

  return;
}
// ============================================================================
} // namespace parameter_matrix
} // namespace nosh
