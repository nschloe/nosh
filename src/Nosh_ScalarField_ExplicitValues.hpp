// @HEADER
//
//    Query routines for a vector potential with explicitly given values.
//    Copyright (C) 2011, 2012  Nico Schl\"omer
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
#ifndef NOSH_SCALARFIELD_EXPLICITVALUES_H_
#define NOSH_SCALARFIELD_EXPLICITVALUES_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Epetra_MultiVector.h>
#include <Teuchos_Array.hpp>

#include "Nosh_ScalarField_Virtual.hpp"
#include "Nosh_StkMesh.hpp"

namespace Nosh {
namespace ScalarField {
class ExplicitValues : public Virtual
{
public:
ExplicitValues(const Nosh::StkMesh &mesh,
               const std::string &fieldName
               );

virtual
~ExplicitValues();

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

//! Get parameter names and initial values.
virtual
const std::map<std::string,double>
getParameters() const;

protected:
private:

Teuchos::ArrayRCP<double> nodeValues_;
};
} // namespace ScalarField
} // namespace Nosh
#endif // NOSH_SCALARFIELD_EXPLICITVALUES_H_
