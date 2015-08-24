// @HEADER
//
//    Query routines for a vector field.
//    Copyright (C) 2012  Nico Schlömer
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
#ifndef NOSH_VECTORFIELD_BASE_H_
#define NOSH_VECTORFIELD_BASE_H_

#include <string>
#include <map>

#include <Eigen/Dense>

namespace nosh
{
namespace vector_field
{
class base
{
public:
  virtual
  ~base(){};

  virtual
  void
  set_parameters(const std::map<std::string, double> & params) = 0;

  //! get parameter names and initial values.
  virtual
  const std::map<std::string, double>
  get_parameters() const = 0;

  //! Projection of the vector field onto an edge at the midpoint of the edge.
  virtual
  double
  get_edge_projection(const unsigned int edge_index) const = 0;

  virtual
  double
  get_d_edge_projection_dp(
      const unsigned int edge_index,
      const std::string & dParamName
      ) const = 0;

  virtual
  Eigen::Vector3d
  eval(const Eigen::Vector3d & x) const = 0;

  virtual
  unsigned int
  degree() const = 0;

protected:
private:
};
} // namespace vector_field
} // namespace nosh
#endif // NOSH_VECTORFIELD_BASE_H_
