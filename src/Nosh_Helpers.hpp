// @HEADER
//
//    Helpers functions.
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
#ifndef NOSH_HELPERS_H
#define NOSH_HELPERS_H
// =============================================================================
// includes
#include <string>

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>
// =============================================================================
namespace Nosh {
namespace Helpers {
void
StkMeshRead(const Epetra_Comm &comm,
            const std::string &fileName,
            const int step,
            Teuchos::ParameterList &data
            );
} // namespace Helpers
} // namespace Nosh
// =============================================================================
#endif // NOSH_HELPERS_H
