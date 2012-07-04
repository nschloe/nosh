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
#ifndef NOSH_SCALARPOTENTIAL_CONSTANT_H_
#define NOSH_SCALARPOTENTIAL_CONSTANT_H_
// =============================================================================
// forward defs
class Epetra_Vector;
class Epetra_Map;
// =============================================================================
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Epetra_Comm.h>

#include "Nosh_ScalarPotential_Virtual.hpp"
// =============================================================================
namespace Nosh {
namespace ScalarPotential {
class Constant: public Virtual
{
public:
Constant(const double alpha
         );

Epetra_Vector
createPInit_(const Epetra_Map & map);

~Constant();

virtual
double
getV(const unsigned int nodeIndex,
     const Teuchos::Array<double> & p
     ) const;

virtual
double
getdVdP(const unsigned int nodeIndex,
        const unsigned int parameterIndex,
        const Teuchos::Array<double> & p
        ) const;

//! Gets the parameters from this module.
virtual
Teuchos::RCP<const Teuchos::Array<double> >
get_p_init() const;

//! Get the parameter names.
virtual
Teuchos::RCP<const Teuchos::Array<std::string> >
get_p_names() const;

protected:
private:

const double alpha_;

};
} // namespace ScalarPotential
} // namespace Nosh
#endif // NOSH_SCALARPOTENTIAL_CONSTANT_H_
