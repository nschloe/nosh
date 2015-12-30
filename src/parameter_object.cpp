// @HEADER
//
//    parameter_object class for matrix constructors.
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

#include "parameter_object.hpp"

namespace nosh
{
// ============================================================================
parameter_object::
parameter_object():
  build_parameters_()
{
}
// ============================================================================
parameter_object::
~parameter_object()
{
}
// ============================================================================
void
parameter_object::
set_parameters(const std::map<std::string, double> &params)
{
  std::cout << ">> po::set_parameters" << std::endl;
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

  std::cout << "aa " << needs_refill << std::endl;

  if (needs_refill) {
    std::cout << "about to refill" << std::endl;
    this->refill_(params);
    std::cout << "just refilled" << std::endl;
    for (auto &build_param: build_parameters_) {
      build_param.second = params.at(build_param.first);
    }
  }

  std::cout << "   po::set_parameters >>" << std::endl;
  return;
}
// ============================================================================
} // namespace nosh
