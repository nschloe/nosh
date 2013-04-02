// @HEADER
//
//    Nosh virtual model evaluator.
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

#include "nosh/ModelEvaluator_Virtual.hpp"

namespace Nosh {
namespace ModelEvaluator {
// ============================================================================
Virtual::
Virtual ()
{
}
// ============================================================================
Virtual::
~Virtual()
{
}
// =============================================================================
double
Virtual::
norm(const Epetra_Vector &psi) const
{
  return sqrt(this->innerProduct(psi, psi));
}
// =============================================================================
} // namespace ModelEvaluator
} // namespace Nosh