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
#include <map>
#include <string>

#include <Teuchos_RCP.hpp>

#include "ScalarField_Virtual.hpp"

#include "Mesh.hpp"
// =============================================================================
namespace Nosh
{
namespace ScalarField
{
class Constant: public Virtual
{
public:
  Constant(
      const Nosh::Mesh & mesh,
      const double c,
      const std::string & param1Name = "",
      const double param1InitValue = 0.0
      );

  Tpetra::Vector<double,int,int>
  createPInit_(const Tpetra::Map<int,int> & map);

  ~Constant();

  //! Get the parameter names and intial values.
  virtual
  const std::map<std::string, double>
  getParameters() const;

  virtual
  const Tpetra::Vector<double,int,int>
  getV(const std::map<std::string, double> & params
     ) const;

  virtual
  const Tpetra::Vector<double,int,int>
  getdVdP(const std::map<std::string, double> & params,
          const std::string & paramName
        ) const;

protected:
private:
  const std::shared_ptr<const Tpetra::Map<int,int>> map_;
  const double c_;
  const std::string param1Name_;
  const double param1InitValue_;
};
} // namespace ScalarField
} // namespace Nosh
#endif // NOSH_SCALARFIELD_CONSTANT_H_
