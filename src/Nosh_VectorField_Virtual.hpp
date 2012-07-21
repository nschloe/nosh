// @HEADER
//
//    Query routines for a vector field.
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
#ifndef NOSH_VECTORFIELD_VIRTUAL_H_
#define NOSH_VECTORFIELD_VIRTUAL_H_
// =============================================================================
// forward decls
class Epetra_Vector;
class Epetra_Map;
// =============================================================================
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
// =============================================================================
namespace Nosh {
namespace VectorField {
class Virtual
{
public:
Virtual();

virtual
~Virtual();

//! Projection of the vector field onto an edge at the midpoint of the edge.
virtual
double
getEdgeProjection(const unsigned int edgeIndex,
                  const Teuchos::Array<double> &params
                  ) const = 0;

virtual
double
getDEdgeProjectionDp(const unsigned int edgeIndex,
                     const Teuchos::Array<double> &params,
                     const unsigned int parameterIndex
                     ) const = 0;

//! Gets the current parameters from this module.
virtual
Teuchos::RCP<const Teuchos::Array<double> >
get_p_init() const = 0;

//! Get the parameter names.
virtual
Teuchos::RCP<const Teuchos::Array<std::string> >
get_p_names() const = 0;

protected:
private:
};
} // namespace VectorField
} // namespace Nosh
#endif // NOSH_VECTORFIELD_VIRTUAL_H_
