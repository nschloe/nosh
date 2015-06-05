// @HEADER
//
//    Query routines for a virtual scalar field.
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
#ifndef NOSH_SCALARFIELD_VIRTUAL_H_
#define NOSH_SCALARFIELD_VIRTUAL_H_
// =============================================================================
#include <map>
#include <string>

#include <Tpetra_Vector.hpp>

namespace Nosh
{
namespace ScalarField
{
class Virtual
{
public:
  Virtual();

  virtual
  ~Virtual();

  virtual
  const Tpetra::Vector<double,int,int>
  getV(const std::map<std::string, double> & params
     ) const = 0;

  virtual
  const Tpetra::Vector<double,int,int>
  getdVdP(const std::map<std::string, double> & params,
          const std::string & paramName
        ) const = 0;

  //! Get parameter names and initial values.
  virtual
  const std::map<std::string, double>
  getParameters() const = 0;

protected:
private:
};
} // namespace ScalarField
} // namespace Nosh
#endif // NOSH_SCALARFIELD_VIRTUAL_H_
