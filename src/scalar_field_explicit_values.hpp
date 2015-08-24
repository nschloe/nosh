// @HEADER
//
//    Query routines for a vector potential with explicitly given values.
//    Copyright (C) 2011, 2012  Nico Schlömer
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

#include <map>
#include <string>

#include <Teuchos_RCP.hpp>
#include <Tpetra_Vector.hpp>

#include "scalar_field_base.hpp"
#include "mesh.hpp"

namespace nosh
{
namespace scalar_field
{
class explicit_values : public base
{
public:
  explicit_values(
      const nosh::mesh &mesh,
      const std::string &field_name
      );

  virtual
  ~explicit_values();

//! Get parameter names and initial values.
  virtual
  const std::map<std::string, double>
  get_parameters() const;

  virtual
  const Tpetra::Vector<double,int,int>
  get_v(const std::map<std::string, double> & params) const;

  virtual
  const Tpetra::Vector<double,int,int>
  get_dvdp(
      const std::map<std::string, double> & params,
      const std::string & param_name
      ) const;

protected:
private:
  const std::shared_ptr<const Tpetra::Vector<double,int,int>> node_values_;
};
} // namespace scalar_field
} // namespace nosh
#endif // NOSH_SCALARFIELD_EXPLICITVALUES_H_
