// @HEADER
//
//    Custom scalar potential.
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
#ifndef MYSCALARFIELD_H_
#define MYSCALARFIELD_H_

// forward defs
class Epetra_Vector;
class Epetra_Map;
namespace Nosh{
class StkMesh;
}

#include <map>
#include <memory>
#include <string>

#include <Epetra_Comm.h>

#include "nosh/ScalarField_Virtual.hpp"

class MyScalarField: public Nosh::ScalarField::Virtual
{
public:
MyScalarField(const std::shared_ptr<const Nosh::StkMesh> & mesh);

Epetra_Vector
createPInit_(const Epetra_Map & map);

~MyScalarField();

//! Get parameter names and initial values.
virtual
const std::map<std::string, double>
getParameters() const;

virtual
const Epetra_Vector
getV(const std::map<std::string, double> & params) const;

virtual
const Epetra_Vector
getdVdP(const std::map<std::string, double> & params,
        const std::string & paramName
        ) const;

protected:
private:
const std::shared_ptr<const Nosh::StkMesh> mesh_;
};
#endif // MYSCALARFIELD_H_
