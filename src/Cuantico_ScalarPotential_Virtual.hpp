// @HEADER
//
//    Query routines for the magnetic vector potential.
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
#ifndef CUANTICO_SCALARPOTENTIAL_VIRTUAL_H_
#define CUANTICO_SCALARPOTENTIAL_VIRTUAL_H_
// =============================================================================
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
// =============================================================================
namespace Cuantico {
namespace ScalarPotential {
class Virtual
{
public:
Virtual();

~Virtual();

virtual
double
getV(const unsigned int nodeIndex,
     const Teuchos::Array<double> & p
    ) const = 0;

virtual
double
getdVdP(const unsigned int nodeIndex,
        const unsigned int parameterIndex,
        const Teuchos::Array<double> & p
        ) const = 0;

//! Gets the parameters from this module.
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
} // namespace ScalarPotential
} // namespace Cuantico
#endif // CUANTICO_SCALARPOTENTIAL_VIRTUAL_H_
