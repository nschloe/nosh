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
#ifndef NOSH_SCALARFIELD_CONSTANT_H_
#define NOSH_SCALARFIELD_CONSTANT_H_
// =============================================================================
// forward defs
class Epetra_Vector;
class Epetra_Map;
// =============================================================================
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Epetra_Comm.h>

#include "Nosh_ScalarField_Virtual.hpp"
// =============================================================================
namespace Nosh {
namespace ScalarField {
class Constant: public Virtual
{
public:
Constant(const double c,
         const std::string & param1Name = "",
         const double param1InitValue = 0.0
         );

Epetra_Vector
createPInit_(const Epetra_Map & map);

~Constant();

//! Get the parameter names and intial values.
virtual
const std::map<std::string,double>
getParameters() const;

virtual
double
getV(const unsigned int nodeIndex,
     const std::map<std::string,double> & params
     ) const;

virtual
double
getdVdP(const unsigned int nodeIndex,
        const std::map<std::string,double> & params,
        const std::string & paramName
        ) const;

protected:
private:

const double c_;
const std::string param1Name_;
const double param1InitValue_;

};
} // namespace ScalarField
} // namespace Nosh
#endif // NOSH_SCALARFIELD_CONSTANT_H_
